"""Schema Updater for SpecificTermExtractionSchema in llm_curation_schema.py.

Safely locates SpecificTermExtractionSchema field definitions and updates Literal[...]
type annotations with newly approved ontology candidate terms.
"""

import ast
import contextlib
import difflib
import hashlib
import os
import re
import tempfile
from pathlib import Path
from typing import Literal, get_args, get_origin

from curation_tools.llm_curation.curation_run_store import (
    CurationRunStore,
    SchemaOperation,
)

DEFAULT_SCHEMA_PATH = Path(__file__).resolve().parent / "llm_curation_schema.py"


class SchemaUpdateConflictError(RuntimeError):
    """Raised when the schema changed after an update was rendered."""


@contextlib.contextmanager
def _schema_update_lock(schema_file_path: Path):
    """Serialize cooperating schema writers with an advisory process lock."""
    try:
        import fcntl
    except ImportError as exc:  # pragma: no cover - unsupported on Windows
        raise RuntimeError(
            "Schema updates require a platform file-lock implementation"
        ) from exc

    lock_path = schema_file_path.with_name(f".{schema_file_path.name}.lock")
    with lock_path.open("a+", encoding="utf-8") as lock_file:
        fcntl.flock(lock_file.fileno(), fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(lock_file.fileno(), fcntl.LOCK_UN)


def get_existing_literals(
    schema_file_path: Path = DEFAULT_SCHEMA_PATH,
) -> dict[str, list[str]]:
    """Return a mapping of field_name -> list of string values for Literal fields in SpecificTermExtractionSchema."""
    from curation_tools.llm_curation.llm_curation_schema import (
        SpecificTermExtractionSchema,
    )

    result: dict[str, list[str]] = {}
    for field_name, field_info in SpecificTermExtractionSchema.model_fields.items():
        annotation = field_info.annotation

        def extract_literals(tp):
            origin = get_origin(tp)
            if origin is Literal:
                return list(get_args(tp))
            args = get_args(tp)
            res = []
            for arg in args:
                if get_origin(arg) is Literal:
                    res.extend(get_args(arg))
            return res

        literal_vals = extract_literals(annotation)
        if literal_vals:
            result[field_name] = [str(v) for v in literal_vals if isinstance(v, str)]
    return result


def update_literal_string_in_code(
    lines: list[str], start_line: int, end_line: int, new_terms: list[str]
) -> list[str]:
    """Given lines of code corresponding to an AnnAssign node, find Literal[...] and insert new_terms."""
    segment = "".join(lines[start_line - 1 : end_line])

    # Match Literal[...] spanning single or multiple lines
    match = re.search(r"Literal\s*\[(.*?)\]", segment, re.DOTALL)
    if not match:
        return lines

    literal_inner = match.group(1)

    # Extract existing literal items (strings inside quotes)
    # Using regex to find all quoted strings
    existing_quotes = re.findall(
        r'(?:"([^"\\]*(?:\\.[^"\\]*)*)"|\'([^\'\\]*(?:\\.[^\'\\]*)*)\')', literal_inner
    )
    existing_items = [q[0] if q[0] != "" else q[1] for q in existing_quotes]

    # Filter out terms that are already present (case-exact)
    terms_to_add = [t for t in new_terms if t not in existing_items]
    if not terms_to_add:
        return lines

    # Determine if "Other" is present
    has_other = "Other" in existing_items

    if has_other:
        # Insert before "Other"
        idx = existing_items.index("Other")
        updated_items = existing_items[:idx] + terms_to_add + existing_items[idx:]
    else:
        updated_items = existing_items + terms_to_add

    # Check formatting of the original Literal block (multi-line vs single-line)
    if "\n" in literal_inner:
        # Detect indentation from the lines in literal_inner
        inner_lines = [l for l in literal_inner.split("\n") if l.strip()]
        indent = "            "  # default 12 spaces
        if inner_lines:
            first_line = inner_lines[0]
            indent_match = re.match(r"^(\s*)", first_line)
            if indent_match and indent_match.group(1):
                indent = indent_match.group(1)

        formatted_items = ",\n".join(f'{indent}"{item}"' for item in updated_items)
        new_literal = (
            f"Literal[\n{formatted_items},\n{indent[:-4] if len(indent) >= 4 else ''}]"
        )
    else:
        formatted_items = ", ".join(f'"{item}"' for item in updated_items)
        new_literal = f"Literal[{formatted_items}]"

    new_segment = segment[: match.start()] + new_literal + segment[match.end() :]

    new_lines = (
        lines[: start_line - 1]
        + new_segment.splitlines(keepends=True)
        + lines[end_line:]
    )
    return new_lines


def _generate_updated_schema_code_from_source(
    approved_terms: dict[str, list[str]], schema_file_path: Path, code: str
) -> str:
    """Render an updated schema from one already-read source snapshot."""
    tree = ast.parse(code)

    # Find ClassDef for SpecificTermExtractionSchema
    target_class = None
    for node in tree.body:
        if (
            isinstance(node, ast.ClassDef)
            and node.name == "SpecificTermExtractionSchema"
        ):
            target_class = node
            break

    if not target_class:
        raise ValueError(
            f"Class 'SpecificTermExtractionSchema' not found in {schema_file_path}"
        )

    # Build field_name -> AST node mapping for fields in SpecificTermExtractionSchema
    field_nodes: dict[str, ast.AnnAssign] = {}
    for stmt in target_class.body:
        if isinstance(stmt, ast.AnnAssign) and isinstance(stmt.target, ast.Name):
            field_nodes[stmt.target.id] = stmt

    lines = code.splitlines(keepends=True)

    # Sort field updates in reverse order of line numbers so line replacements don't shift earlier line numbers
    field_updates = []
    for field_name, new_terms in approved_terms.items():
        if not new_terms or field_name not in field_nodes:
            continue
        node = field_nodes[field_name]
        field_updates.append((node.lineno, node.end_lineno, field_name, new_terms))

    field_updates.sort(key=lambda x: x[0], reverse=True)

    for start_line, end_line, field_name, new_terms in field_updates:
        lines = update_literal_string_in_code(lines, start_line, end_line, new_terms)

    return "".join(lines)


def get_schema_diff(
    approved_terms: dict[str, list[str]], schema_file_path: Path = DEFAULT_SCHEMA_PATH
) -> str:
    """Return unified diff string between current schema and updated schema."""
    original_code = schema_file_path.read_text(encoding="utf-8")
    updated_code = _generate_updated_schema_code_from_source(
        approved_terms, schema_file_path, original_code
    )

    diff = difflib.unified_diff(
        original_code.splitlines(keepends=True),
        updated_code.splitlines(keepends=True),
        fromfile=f"a/{schema_file_path.name}",
        tofile=f"b/{schema_file_path.name}",
    )
    return "".join(diff)


def apply_schema_update(
    approved_terms: dict[str, list[str]],
    schema_path: Path = DEFAULT_SCHEMA_PATH,
    run_store: CurationRunStore | None = None,
) -> SchemaOperation:
    """Journal and atomically apply approved schema terms for one native run."""
    schema_path = Path(schema_path).resolve()
    if run_store is None:
        raise ValueError(
            "Schema updates require a CurationRunStore for durable journaling"
        )

    with _schema_update_lock(schema_path):
        original_code = schema_path.read_text(encoding="utf-8")
        updated_code = _generate_updated_schema_code_from_source(
            approved_terms, schema_path, original_code
        )
        diff = "".join(
            difflib.unified_diff(
                original_code.splitlines(keepends=True),
                updated_code.splitlines(keepends=True),
                fromfile=f"a/{schema_path.name}",
                tofile=f"b/{schema_path.name}",
            )
        )
        old_hash = _text_hash(original_code)
        operation = run_store.begin_schema_operation(
            schema_path,
            original_code,
            updated_code,
            diff,
            approved_terms,
        )

        temporary_path: Path | None = None
        try:
            ast.parse(updated_code, filename=str(schema_path))
            compile(updated_code, str(schema_path), "exec")
            fd, temp_name = tempfile.mkstemp(
                prefix=f".{schema_path.name}.",
                suffix=".tmp",
                dir=schema_path.parent,
            )
            temporary_path = Path(temp_name)
            with os.fdopen(fd, "w", encoding="utf-8") as temporary:
                temporary.write(updated_code)
                temporary.flush()
                os.fsync(temporary.fileno())

            if _file_hash(schema_path) != old_hash:
                raise SchemaUpdateConflictError(
                    "Schema changed after rendering; update was not applied"
                )
            os.replace(temporary_path, schema_path)
            temporary_path = None
            _fsync_directory(schema_path.parent)
            run_store.finish_schema_operation(operation.operation_id, "applied")
            operation = SchemaOperation(
                operation.operation_id,
                "applied",
                operation.old_hash,
                operation.new_hash,
            )
        except Exception as exc:
            run_store.finish_schema_operation(
                operation.operation_id, "failed", str(exc)
            )
            raise
        finally:
            if temporary_path is not None:
                temporary_path.unlink(missing_ok=True)

    return operation


def _file_hash(file_path: Path) -> str | None:
    """Return a SHA-256 hash for a file when it exists."""
    if not file_path.is_file():
        return None
    return hashlib.sha256(file_path.read_bytes()).hexdigest()


def _text_hash(text: str) -> str:
    """Return a SHA-256 hash for rendered schema text."""
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _fsync_directory(directory: Path) -> None:
    """Flush a directory entry when the platform supports directory fsync."""
    try:
        directory_fd = os.open(directory, os.O_RDONLY)
    except OSError:
        return
    try:
        os.fsync(directory_fd)
    finally:
        os.close(directory_fd)
