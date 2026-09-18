import json
import logging
import argparse
import math
import os
from pathlib import Path

import torch
from datasets import Dataset
from peft import LoraConfig, get_peft_model, TaskType, prepare_model_for_kbit_training
from transformers import (
    AutoModelForCausalLM,
    AutoTokenizer,
    BitsAndBytesConfig,
    TrainingArguments,
)
from trl import SFTTrainer

log = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s"
)


PROMPT_TEMPLATE = (
    "### Instruction:\n{instruction}\n\n"
    "### Input:\n{input}\n\n"
    "### Response:\n{output}"
)

# At inference time, the model is given everything up to and including
INFERENCE_TEMPLATE = (
    "### Instruction:\n{instruction}\n\n" "### Input:\n{input}\n\n" "### Response:\n"
)


def format_for_sft(record):
    """
    Format a training record as a single string for SFTTrainer.

    The metadata field is intentionally excluded — it contains dataset
    labels (disease, modality, fitness_class) that must not be seen
    by the model during training. The model should learn to predict
    the output from instruction + input alone.

    Parameters
    ----------
    record : dict
        Training record with instruction, input, output fields.

    Returns
    -------
    str — full prompt including the target response.
    """
    return PROMPT_TEMPLATE.format(
        instruction=record["instruction"],
        input=record["input"],
        output=record["output"],
    )


def load_split(jsonl_path):
    """
    Load records from a JSONL file and return as a HuggingFace Dataset.

    Validates that each record has instruction, input, output fields.
    Skips malformed records with a warning.

    Parameters
    ----------
    jsonl_path : str or Path

    Returns
    -------
    datasets.Dataset with a single 'text' column containing
    formatted SFT prompts.
    """
    records = []
    skipped = 0

    with open(jsonl_path) as f:
        for i, line in enumerate(f, start=1):
            line = line.strip()
            if not line:
                continue
            try:
                record = json.loads(line)
            except json.JSONDecodeError as e:
                log.warning(f"Line {i}: JSON parse error — {e}. Skipping.")
                skipped += 1
                continue

            required = {"instruction", "input", "output"}
            if not required.issubset(record.keys()):
                log.warning(
                    f"Line {i}: missing fields {required - set(record.keys())}. Skipping."
                )
                skipped += 1
                continue

            records.append({"text": format_for_sft(record)})

    if skipped:
        log.warning(f"Skipped {skipped} malformed records in {jsonl_path}")

    log.info(f"Loaded {len(records)} records from {jsonl_path}")
    return Dataset.from_list(records)


def load_splits(splits_dir):
    """
    Load train and val splits from a directory created by pipeline.py --split.

    Expects:
        splits_dir/train.jsonl
        splits_dir/val.jsonl
        splits_dir/split_manifest.json

    Parameters
    ----------
    splits_dir : str or Path

    Returns
    -------
    tuple of (train_dataset, val_dataset, split_manifest)
    """
    splits_dir = Path(splits_dir)

    train_path = splits_dir / "train.jsonl"
    val_path = splits_dir / "val.jsonl"
    manifest_path = splits_dir / "split_manifest.json"

    for path in [train_path, val_path, manifest_path]:
        if not path.exists():
            raise FileNotFoundError(
                f"Expected split file not found: {path}\n"
                f"Run: python pipeline.py --split --corpus <corpus.jsonl> "
                f"--output_dir {splits_dir}"
            )

    train_dataset = load_split(train_path)
    val_dataset = load_split(val_path)

    with open(manifest_path) as f:
        split_manifest = json.load(f)

    log.info(
        f"Splits loaded: {len(train_dataset)} train / {len(val_dataset)} val records"
    )
    log.info(f"Split manifest: {len(split_manifest)} genes assigned to splits")

    return train_dataset, val_dataset, split_manifest


def build_lora_model(
    model_name,
    use_qlora=False,
    lora_r=16,
    lora_alpha=32,
    lora_dropout=0.05,
    max_seq_length=2048,
):
    """
    Load a causal LM and wrap it with LoRA (or QLoRA) adapters.

    LoRA freezes the base model weights and adds small trainable
    rank-decomposition matrices to the attention projection layers.
    This reduces trainable parameters from ~2.7B to ~10-20M while
    preserving most of the base model's knowledge.

    QLoRA additionally quantises the frozen base model to 4-bit,
    reducing GPU memory from ~16GB to ~6GB for BioMedLM 2.7B.

    Parameters
    ----------
    model_name : str
        HuggingFace model ID, e.g. "stanford-crfm/BioMedLM".
    use_qlora : bool
        If True, load base model in 4-bit NF4 quantisation (QLoRA).
    lora_r : int
        LoRA rank. Higher = more parameters, more capacity.
    lora_alpha : int
        LoRA scaling factor. Convention: alpha = 2 * r.
    lora_dropout : float
        Dropout applied to LoRA layers. 0.05 is a safe default.
    max_seq_length : int
        Maximum token sequence length. 2048 is safe for BioMedLM.

    Returns
    -------
    tuple of (model, tokenizer)
    """
    log.info(f"Loading tokenizer from {model_name}...")
    tokenizer = AutoTokenizer.from_pretrained(model_name)

    # BioMedLM (GPT-2 architecture) has no pad token by default
    # SFTTrainer requires one for batched training
    if tokenizer.pad_token is None:
        tokenizer.pad_token = tokenizer.eos_token
        log.info("Set pad_token = eos_token (required for batched SFT)")

    tokenizer.model_max_length = max_seq_length

    bnb_config = None
    if use_qlora:
        if not torch.cuda.is_available():
            raise RuntimeError(
                "QLoRA requires a CUDA GPU. "
                "Run without --use_qlora on CPU or use a GPU environment."
            )
        bnb_config = BitsAndBytesConfig(
            load_in_4bit=True,
            bnb_4bit_quant_type="nf4",
            bnb_4bit_compute_dtype=torch.bfloat16,
            bnb_4bit_use_double_quant=True,
        )
        log.info("QLoRA enabled: loading base model in 4-bit NF4")

    log.info(f"Loading model from {model_name}...")
    model = AutoModelForCausalLM.from_pretrained(
        model_name,
        quantization_config=bnb_config,
        device_map="auto" if torch.cuda.is_available() else None,
        torch_dtype=torch.bfloat16 if torch.cuda.is_available() else torch.float32,
        trust_remote_code=True,  # required for BioMedLM
    )

    if use_qlora:
        model = prepare_model_for_kbit_training(model)

    # BioMedLM is GPT-2 architecture — uses c_attn/c_proj not q_proj/k_proj
    if any(name in model_name.lower() for name in ["gpt2", "gpt-2", "biomedlm"]):
        target_modules = ["c_attn", "c_proj"]
    else:
        target_modules = ["q_proj", "k_proj", "v_proj", "o_proj"]

    # q_proj, k_proj, v_proj, o_proj are the query/key/value/output
    lora_config = LoraConfig(
        task_type=TaskType.CAUSAL_LM,
        r=lora_r,
        lora_alpha=lora_alpha,
        lora_dropout=lora_dropout,
        target_modules=target_modules,
        bias="none",
    )

    model = get_peft_model(model, lora_config)

    trainable_params = sum(p.numel() for p in model.parameters() if p.requires_grad)
    total_params = sum(p.numel() for p in model.parameters())
    log.info(
        f"Trainable parameters: {trainable_params:,} / {total_params:,} "
        f"({100 * trainable_params / total_params:.2f}%)"
    )

    return model, tokenizer


def build_training_args(
    output_dir,
    epochs=3,
    batch_size=4,
    grad_accum=4,
    lr=2e-4,
    wandb_project=None,
):
    """
    Build HuggingFace TrainingArguments for SFT.

    Effective batch size = batch_size * grad_accum = 4 * 4 = 16.

    Cosine learning rate schedule: starts at lr, decays to ~0 by end
    of training. Better than linear for fine-tuning — avoids sharp
    drops early in training.

    Parameters
    ----------
    output_dir : str
        Directory to save checkpoints and final adapter.
    epochs : int
        Number of training epochs.
    batch_size : int
        Per-device batch size.
    grad_accum : int
        Gradient accumulation steps.
        Effective batch = batch_size * grad_accum = 16.
    lr : float
        Peak learning rate.
    wandb_project : str or None
        W&B project name. If None and WANDB_API_KEY is not set,
        W&B logging is disabled gracefully.

    Returns
    -------
    transformers.TrainingArguments
    """

    use_wandb = (
        wandb_project is not None and os.environ.get("WANDB_API_KEY") is not None
    )
    report_to = "wandb" if use_wandb else "none"

    if wandb_project and not use_wandb:
        log.warning(
            "WANDB_API_KEY not set — W&B logging disabled. "
            "Set the environment variable to enable experiment tracking."
        )

    return TrainingArguments(
        output_dir=output_dir,
        num_train_epochs=epochs,
        per_device_train_batch_size=batch_size,
        per_device_eval_batch_size=batch_size,
        gradient_accumulation_steps=grad_accum,
        learning_rate=lr,
        lr_scheduler_type="cosine",
        warmup_ratio=0.05,
        weight_decay=0.01,
        fp16=False,
        bf16=torch.cuda.is_available(),
        logging_steps=10,
        eval_strategy="epoch",
        save_strategy="epoch",
        load_best_model_at_end=True,
        metric_for_best_model="eval_loss",
        greater_is_better=False,
        report_to=report_to,
        run_name=wandb_project,
        dataloader_num_workers=0,
        remove_unused_columns=False,
    )


def compute_perplexity(eval_loss):
    """
    Compute perplexity from cross-entropy eval loss.

    Perplexity = exp(loss). Lower is better.
    A perplexity of 1.0 means the model perfectly predicts every token.
    A perplexity of 100 means the model is very uncertain.

    This is the primary convergence metric logged to W&B per epoch.

    Parameters
    ----------
    eval_loss : float

    Returns
    -------
    float
    """
    return math.exp(eval_loss)


def train(
    model,
    tokenizer,
    train_dataset,
    val_dataset,
    training_args,
    split_manifest,
    max_seq_length=2048,
):
    """
    Run supervised fine-tuning with SFTTrainer.

    SFTTrainer (from TRL) handles:
    - Packing short sequences together for efficiency
    - Masking the prompt tokens so loss is only computed on the response
    - Gradient checkpointing for memory efficiency

    After training, saves:
    - LoRA adapter weights to output_dir/final_adapter/
    - Split manifest to output_dir/split_manifest.json
      (so evaluate.py can reconstruct the exact same test genes)

    Parameters
    ----------
    model : peft.PeftModel
    tokenizer : transformers.PreTrainedTokenizer
    train_dataset : datasets.Dataset
    val_dataset : datasets.Dataset
    training_args : transformers.TrainingArguments
    split_manifest : dict — gene → split label
    max_seq_length : int

    Returns
    -------
    trainer : trl.SFTTrainer (contains training history)
    """

    trainer = SFTTrainer(
        model=model,
        processing_class=tokenizer,
        train_dataset=train_dataset,
        eval_dataset=val_dataset,
        args=training_args,
    )

    log.info("Starting training...")
    log.info(
        f"  Epochs: {training_args.num_train_epochs} | "
        f"  Batch: {training_args.per_device_train_batch_size} × "
        f"{training_args.gradient_accumulation_steps} = "
        f"{training_args.per_device_train_batch_size * training_args.gradient_accumulation_steps} effective | "
        f"  LR: {training_args.learning_rate}"
    )

    trainer.train()

    eval_results = trainer.evaluate()
    final_perplexity = compute_perplexity(eval_results["eval_loss"])
    log.info(f"Final val loss: {eval_results['eval_loss']:.4f}")
    log.info(f"Final val perplexity: {final_perplexity:.2f}")

    if training_args.report_to == "wandb":
        try:
            import wandb

            wandb.log({"final_val_perplexity": final_perplexity})
        except ImportError:
            pass

    adapter_path = Path(training_args.output_dir) / "final_adapter"
    adapter_path.mkdir(parents=True, exist_ok=True)
    model.save_pretrained(str(adapter_path))
    tokenizer.save_pretrained(str(adapter_path))
    log.info(f"Saved LoRA adapter to {adapter_path}")

    manifest_path = Path(training_args.output_dir) / "split_manifest.json"
    with open(manifest_path, "w") as f:
        json.dump(split_manifest, f, indent=2)
    log.info(f"Saved split manifest to {manifest_path}")

    return trainer


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Fine-tune BioMedLM on perturbation training records using LoRA. "
            "Run prepare_training_data.py and pipeline.py --split first."
        )
    )

    # Data
    parser.add_argument(
        "--splits_dir",
        type=str,
        required=True,
        help=(
            "Directory containing train.jsonl, val.jsonl, split_manifest.json. "
            "Created by: python pipeline.py --split --corpus <file> --output_dir <dir>"
        ),
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        required=True,
        help="Directory to save checkpoints and final LoRA adapter",
    )

    # Model
    parser.add_argument(
        "--model_name",
        type=str,
        default="stanford-crfm/BioMedLM",
        help="HuggingFace model ID (default: stanford-crfm/BioMedLM)",
    )
    parser.add_argument(
        "--use_qlora",
        action="store_true",
        help="Use QLoRA (4-bit quantisation) — reduces GPU memory at slight quality cost",
    )
    parser.add_argument(
        "--max_seq_length",
        type=int,
        default=2048,
        help="Maximum token sequence length (default: 2048)",
    )

    # LoRA hyperparameters
    parser.add_argument(
        "--lora_r",
        type=int,
        default=16,
        help="LoRA rank (default: 16)",
    )
    parser.add_argument(
        "--lora_alpha",
        type=int,
        default=32,
        help="LoRA alpha scaling factor (default: 32 = 2 × lora_r)",
    )
    parser.add_argument(
        "--lora_dropout",
        type=float,
        default=0.05,
        help="LoRA dropout (default: 0.05)",
    )

    # Training hyperparameters
    parser.add_argument(
        "--epochs",
        type=int,
        default=3,
        help="Number of training epochs (default: 3)",
    )
    parser.add_argument(
        "--batch_size",
        type=int,
        default=4,
        help="Per-device batch size (default: 4)",
    )
    parser.add_argument(
        "--grad_accum",
        type=int,
        default=4,
        help="Gradient accumulation steps (default: 4 → effective batch 16)",
    )
    parser.add_argument(
        "--lr",
        type=float,
        default=2e-4,
        help="Peak learning rate (default: 2e-4)",
    )

    # Experiment tracking
    parser.add_argument(
        "--wandb_project",
        type=str,
        default=None,
        help="W&B project name. Requires WANDB_API_KEY env var to be set.",
    )

    args = parser.parse_args()

    train_dataset, val_dataset, split_manifest = load_splits(args.splits_dir)

    model, tokenizer = build_lora_model(
        model_name=args.model_name,
        use_qlora=args.use_qlora,
        lora_r=args.lora_r,
        lora_alpha=args.lora_alpha,
        lora_dropout=args.lora_dropout,
        max_seq_length=args.max_seq_length,
    )

    training_args = build_training_args(
        output_dir=args.output_dir,
        epochs=args.epochs,
        batch_size=args.batch_size,
        grad_accum=args.grad_accum,
        lr=args.lr,
        wandb_project=args.wandb_project,
    )

    trainer = train(
        model=model,
        tokenizer=tokenizer,
        train_dataset=train_dataset,
        val_dataset=val_dataset,
        training_args=training_args,
        split_manifest=split_manifest,
        max_seq_length=args.max_seq_length,
    )

    print(f"\nTraining complete.")
    print(f"LoRA adapter saved to: {args.output_dir}/final_adapter/")
    print(f"Split manifest saved to: {args.output_dir}/split_manifest.json")
    print(
        f"To evaluate: python evaluate.py "
        f"--adapter_dir {args.output_dir}/final_adapter/ "
        f"--splits_dir {args.splits_dir}"
    )


if __name__ == "__main__":
    main()
