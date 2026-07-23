"""Schema Updater for SpecificTermExtractionSchema in llm_curation_schema.py.

Safely locates SpecificTermExtractionSchema field definitions and updates Literal[...]
type annotations with newly approved ontology candidate terms.
"""

import ast
import difflib
import re
import shutil
from pathlib import Path
from typing import Literal, get_args, get_origin, get_type_hints

DEFAULT_SCHEMA_PATH = Path(__file__).resolve().parent / "llm_curation_schema.py"


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


def generate_updated_schema_code(
    approved_terms: dict[str, list[str]], schema_file_path: Path = DEFAULT_SCHEMA_PATH
) -> str:
    """Generate the updated Python code for llm_curation_schema.py without writing to disk."""
    code = schema_file_path.read_text(encoding="utf-8")
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
    updated_code = generate_updated_schema_code(approved_terms, schema_file_path)

    diff = difflib.unified_diff(
        original_code.splitlines(keepends=True),
        updated_code.splitlines(keepends=True),
        fromfile=f"a/{schema_file_path.name}",
        tofile=f"b/{schema_file_path.name}",
    )
    return "".join(diff)


def apply_schema_update(
    approved_terms: dict[str, list[str]],
    schema_file_path: Path = DEFAULT_SCHEMA_PATH,
    create_backup: bool = True,
) -> str:
    """Apply approved terms to SpecificTermExtractionSchema in llm_curation_schema.py.

    Creates backup if requested, writes changes, and returns diff string.
    """
    if create_backup:
        backup_path = schema_file_path.with_suffix(".py.bak")
        shutil.copy(schema_file_path, backup_path)

    updated_code = generate_updated_schema_code(approved_terms, schema_file_path)
    diff = get_schema_diff(approved_terms, schema_file_path)

    schema_file_path.write_text(updated_code, encoding="utf-8")
    return diff
