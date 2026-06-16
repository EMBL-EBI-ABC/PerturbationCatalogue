import argparse
import html
import re
import xml.etree.ElementTree as ET
from pathlib import Path

try:
    from lxml import etree as LXML_ET
except ImportError:
    LXML_ET = None


PARAGRAPH_TAGS = {"para", "simple-para", "p", "text"}
SECTION_TAGS = {"section", "sections", "sec"}
SECTION_TITLE_TAGS = {"section-title", "title", "title-text"}
LIST_TAGS = {"list", "list-item", "item", "li"}
CAPTION_TAGS = {"caption", "legend", "caption-para"}
ABSTRACT_TAGS = {"abstract"}
REFERENCE_ENTRY_TAGS = {"bib-reference", "reference", "ref", "ref-info"}
REFERENCE_TEXT_TAGS = {"source-text", "ref-fulltext", "mixed-citation", "citation"}
KEYWORD_TAGS = {"keyword", "subject", "term"}
SKIP_BODY_SECTION_TITLES = {"graphical abstract", "keywords"}

def print_status_block(title: str, *lines: str) -> None:
    separator = "=" * 80
    print(f"\n{separator}")
    print(title)
    for line in lines:
        print(line)
    print(separator)

def get_tag_name(elem):
    """Return the local XML tag name without its namespace."""
    return elem.tag.split("}")[-1] if isinstance(elem.tag, str) else ""


def normalize_text(value):
    """Collapse whitespace and decode common HTML/XML entities."""
    if not value:
        return ""
    return re.sub(r"\s+", " ", html.unescape(value)).strip()


def split_keyword_blob(keyword_blob):
    """Split BioC keyword blobs into approximate phrase-level keywords."""
    if not keyword_blob:
        return []

    if re.search(r"[;,|/]", keyword_blob):
        return [normalize_text(part) for part in re.split(r"\s*[;,|/]\s*", keyword_blob) if normalize_text(part)]

    roman_numerals = {"I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X"}
    connector_words = {
        "a",
        "an",
        "and",
        "by",
        "for",
        "in",
        "of",
        "on",
        "or",
        "the",
        "to",
        "via",
        "with",
        "without",
    }

    phrases = []
    current_phrase = []
    for token in keyword_blob.split():
        stripped = token.strip("()[]{}.,;:")
        if not current_phrase:
            current_phrase.append(token)
            continue

        if (
            stripped[:1].islower()
            or stripped.lower() in connector_words
            or stripped in roman_numerals
        ):
            current_phrase.append(token)
            continue

        phrases.append(normalize_text(" ".join(current_phrase)))
        current_phrase = [token]

    if current_phrase:
        phrases.append(normalize_text(" ".join(current_phrase)))

    return [phrase for phrase in phrases if phrase]


def get_node_text(elem):
    """Flatten descendant text into a single normalized string."""
    if elem is None:
        return ""
    return normalize_text("".join(elem.itertext()))


def find_first_ns_agnostic(node, tag_name):
    """Find the first descendant whose local name matches the requested tag."""
    for elem in node.iter():
        if get_tag_name(elem) == tag_name:
            return elem
    return None


def find_all_ns_agnostic(node, tag_name):
    """Find all descendants whose local name matches the requested tag."""
    return [elem for elem in node.iter() if get_tag_name(elem) == tag_name]


def find_first_matching(node, tag_names):
    """Find the first descendant whose local name is in the provided tag set."""
    tag_names = set(tag_names)
    for elem in node.iter():
        if get_tag_name(elem) in tag_names:
            return elem
    return None


def deduplicate_preserve_order(values):
    """Remove duplicates while preserving input order."""
    seen = set()
    deduplicated = []
    for value in values:
        if value and value not in seen:
            deduplicated.append(value)
            seen.add(value)
    return deduplicated


def append_unique(parts, value):
    """Append a non-empty block unless it would repeat the last emitted block."""
    if value and (not parts or parts[-1] != value):
        parts.append(value)


def get_direct_child_text(node, tag_names):
    """Return the first direct-child text for the requested local tag names."""
    tag_names = set(tag_names)
    for child in list(node):
        if get_tag_name(child) in tag_names:
            text = get_node_text(child)
            if text:
                return text, child
    return "", None


def find_elsevier_article_root(root):
    """Return the embedded article subtree for Elsevier XML when available."""
    original_text = find_first_ns_agnostic(root, "originalText")
    if original_text is not None and len(original_text):
        return list(original_text)[0]

    serial_item = find_first_ns_agnostic(root, "serial-item")
    if serial_item is not None:
        return serial_item

    return root


def extract_title(root, article_root=None):
    """Extract a paper title from metadata-first candidates."""
    search_roots = []
    coredata = find_first_ns_agnostic(root, "coredata")
    if coredata is not None:
        search_roots.append(coredata)
    if article_root is not None:
        search_roots.append(article_root)
    search_roots.append(root)

    for search_root in search_roots:
        for tag_names in ({"title"}, {"article-title"}, {"title-text"}):
            candidate = find_first_matching(search_root, tag_names)
            if candidate is not None:
                text = get_node_text(candidate)
                if text:
                    return text

    return "Unknown Title"


def extract_authors(root, article_root=None):
    """Extract author names from creator nodes or structured author elements."""
    authors = []

    for creator in find_all_ns_agnostic(root, "creator"):
        text = get_node_text(creator)
        if text:
            authors.append(text)

    if authors:
        return deduplicate_preserve_order(authors)

    search_roots = [article_root, root] if article_root is not None else [root]
    for search_root in search_roots:
        for author in find_all_ns_agnostic(search_root, "author"):
            given = find_first_ns_agnostic(author, "given-name")
            surname = find_first_ns_agnostic(author, "surname")
            initials = find_first_ns_agnostic(author, "initials")
            parts = [get_node_text(given), get_node_text(surname)]
            name = normalize_text(" ".join(part for part in parts if part))
            if not name:
                name = get_node_text(initials)
            if name:
                authors.append(name)

    return deduplicate_preserve_order(authors)


def extract_keywords(root, article_root=None):
    """Extract keywords from metadata subjects or keyword nodes."""
    keywords = []
    search_roots = [root]
    if article_root is not None:
        search_roots.append(article_root)

    for search_root in search_roots:
        for elem in search_root.iter():
            tag_name = get_tag_name(elem)
            if tag_name in KEYWORD_TAGS:
                text = get_node_text(elem)
                if text:
                    keywords.append(text)

    return deduplicate_preserve_order(keywords)


def render_list(node):
    """Render list-like XML nodes into Markdown bullets."""
    items = []
    for child in list(node):
        tag_name = get_tag_name(child)
        if tag_name in {"list-item", "item", "li"}:
            text = get_node_text(child)
            if text:
                append_unique(items, f"- {text}")
        else:
            items.extend(render_container(child, 0, set(), False))
    return items


def render_container(node, heading_level, skip_section_titles, allow_captions):
    """Render XML content blocks into ordered Markdown paragraphs and headings."""
    parts = []
    for child in list(node):
        tag_name = get_tag_name(child)

        if tag_name in SECTION_TAGS:
            parts.extend(render_section(child, heading_level, skip_section_titles, allow_captions))
            continue

        if tag_name in PARAGRAPH_TAGS:
            text = get_node_text(child)
            if text:
                append_unique(parts, text)
            continue

        if tag_name in LIST_TAGS:
            for item in render_list(child):
                append_unique(parts, item)
            continue

        if allow_captions and tag_name in CAPTION_TAGS:
            text = get_node_text(child)
            if text:
                append_unique(parts, text)
            continue

        parts.extend(render_container(child, heading_level, skip_section_titles, allow_captions))

    return parts


def render_section(node, heading_level, skip_section_titles, allow_captions):
    """Render a semantic section subtree into Markdown with nested headings."""
    parts = []
    title, title_elem = get_direct_child_text(node, SECTION_TITLE_TAGS)
    normalized_title = title.casefold()

    if title and normalized_title not in skip_section_titles:
        append_unique(parts, f"{'#' * min(6, heading_level)} {title}")

    for child in list(node):
        if title_elem is not None and child is title_elem:
            continue

        tag_name = get_tag_name(child)
        if tag_name in SECTION_TAGS:
            parts.extend(render_section(child, heading_level + 1, skip_section_titles, allow_captions))
        elif tag_name in PARAGRAPH_TAGS:
            text = get_node_text(child)
            if text:
                append_unique(parts, text)
        elif tag_name in LIST_TAGS:
            for item in render_list(child):
                append_unique(parts, item)
        elif allow_captions and tag_name in CAPTION_TAGS:
            text = get_node_text(child)
            if text:
                append_unique(parts, text)
        else:
            parts.extend(render_container(child, heading_level + 1, skip_section_titles, allow_captions))

    return parts


def extract_abstract(root, article_root=None):
    """Extract abstract text from semantic abstract nodes or metadata descriptions."""
    search_roots = [article_root, root] if article_root is not None else [root]

    for search_root in search_roots:
        for abstract_node in find_all_ns_agnostic(search_root, "abstract"):
            parts = render_container(abstract_node, 3, set(), False)
            if parts:
                return "\n\n".join(parts)

    for description in find_all_ns_agnostic(root, "description"):
        text = get_node_text(description)
        if text:
            return text

    return ""


def extract_body(root, article_root=None):
    """Extract main body text from article body containers or section trees."""
    search_roots = [article_root, root] if article_root is not None else [root]

    for search_root in search_roots:
        body_node = find_first_ns_agnostic(search_root, "body")
        if body_node is not None:
            parts = render_container(body_node, 2, SKIP_BODY_SECTION_TITLES, True)
            if parts:
                return "\n\n".join(parts)

    for search_root in search_roots:
        for sections_node in find_all_ns_agnostic(search_root, "sections"):
            parts = render_section(sections_node, 2, SKIP_BODY_SECTION_TITLES, True)
            if parts:
                return "\n\n".join(parts)

    return ""


def extract_references(root, article_root=None):
    """Extract reference strings from bibliography-style entries."""
    references = []
    search_roots = [article_root, root] if article_root is not None else [root]

    for search_root in search_roots:
        for entry in search_root.iter():
            if get_tag_name(entry) not in REFERENCE_ENTRY_TAGS:
                continue

            reference_text = ""
            for descendant in entry.iter():
                if get_tag_name(descendant) in REFERENCE_TEXT_TAGS:
                    reference_text = get_node_text(descendant)
                    if reference_text:
                        break

            if not reference_text:
                candidate = get_node_text(entry)
                if len(candidate) > 30:
                    reference_text = candidate

            if reference_text:
                references.append(reference_text)

    return deduplicate_preserve_order(references)


def extract_elsevier_data(root):
    """Extract content from Elsevier XML by prioritizing the embedded article subtree."""
    article_root = find_elsevier_article_root(root)
    return {
        "title": extract_title(root, article_root),
        "authors": extract_authors(root, article_root),
        "keywords": extract_keywords(root, article_root),
        "abstract": extract_abstract(root, article_root),
        "body": extract_body(root, article_root),
        "references": extract_references(root, article_root),
    }


def extract_bioc_data(root):
    """Extract content from BioC XML, such as PubMed Central conversions."""
    title = "Unknown Title"
    authors = []
    keywords = []
    abstract_parts = []
    body_parts = []
    references = []

    for passage in root.findall(".//passage"):
        infons = {infon.get("key"): infon.text for infon in passage.findall("infon")}
        text_elem = passage.find("text")
        text = normalize_text(text_elem.text if text_elem is not None else "")

        sec_type = (infons.get("section_type") or "").upper()
        passage_type = (infons.get("type") or "").lower()

        keyword_blob = normalize_text(infons.get("kwd") or "")
        if keyword_blob:
            keywords.extend(split_keyword_blob(keyword_blob))

        if sec_type == "TITLE":
            if passage_type in {"front", "title"}:
                if text:
                    title = text

            author_items = []
            for key, value in infons.items():
                if key.startswith("name_") and value:
                    parts = dict(part.split(":", 1) for part in value.split(";") if ":" in part)
                    given_names = parts.get("given-names", "")
                    surname = parts.get("surname", "")
                    name = normalize_text(f"{given_names} {surname}")
                    if name:
                        try:
                            order = int(key.split("_", 1)[1])
                        except (IndexError, ValueError):
                            order = len(author_items)
                        author_items.append((order, name))

            authors.extend(name for _, name in sorted(author_items, key=lambda item: item[0]))
            continue

        if not text:
            continue

        elif sec_type == "ABSTRACT":
            if "title" in passage_type:
                if text.casefold() != "abstract":
                    append_unique(abstract_parts, f"### {text}")
            else:
                append_unique(abstract_parts, text)
        elif sec_type == "REF":
            if passage_type == "ref":
                references.append(text)
        elif sec_type == "FIG":
            if passage_type == "fig_caption":
                append_unique(body_parts, f"Figure: {text}")
            else:
                append_unique(body_parts, text)
        elif sec_type == "TABLE":
            if passage_type == "table_caption":
                append_unique(body_parts, f"Table: {text}")
        else:
            if passage_type == "title_1":
                append_unique(body_parts, f"## {text}")
            elif passage_type == "title_2":
                append_unique(body_parts, f"### {text}")
            elif passage_type == "title_3":
                append_unique(body_parts, f"#### {text}")
            else:
                append_unique(body_parts, text)

    return {
        "title": title,
        "authors": deduplicate_preserve_order(authors),
        "keywords": deduplicate_preserve_order(keywords),
        "abstract": "\n\n".join(abstract_parts),
        "body": "\n\n".join(body_parts),
        "references": deduplicate_preserve_order(references),
    }


def extract_regex_fallback(input_file):
    """Fallback extractor for malformed XML when structural parsing fails."""
    content = Path(input_file).read_text(encoding="utf-8", errors="ignore")

    def first_match(patterns):
        for pattern in patterns:
            match = re.search(pattern, content, re.DOTALL | re.IGNORECASE)
            if match:
                text = normalize_text(re.sub(r"<.*?>", " ", match.group(1)))
                if text:
                    return text
        return ""

    title = first_match(
        [
            r"<(?:dc:)?title[^>]*>(.*?)</(?:dc:)?title>",
            r"<article-title[^>]*>(.*?)</article-title>",
            r"<title-text[^>]*>(.*?)</title-text>",
        ]
    ) or "Unknown Title"

    authors = [
        normalize_text(re.sub(r"<.*?>", " ", match.group(1)))
        for match in re.finditer(r"<dc:creator[^>]*>(.*?)</dc:creator>", content, re.DOTALL | re.IGNORECASE)
    ]

    keywords = [
        normalize_text(re.sub(r"<.*?>", " ", match.group(1)))
        for match in re.finditer(r"<dcterms:subject[^>]*>(.*?)</dcterms:subject>", content, re.DOTALL | re.IGNORECASE)
    ]

    abstract = first_match(
        [
            r"<(?:ce:)?abstract[^>]*>(.*?)</(?:ce:)?abstract>",
            r"<dc:description[^>]*>(.*?)</dc:description>",
        ]
    )

    body = first_match([r"<(?:ce:)?body[^>]*>(.*?)</(?:ce:)?body>"])

    references = [
        normalize_text(re.sub(r"<.*?>", " ", match.group(1)))
        for match in re.finditer(
            r"<(?:ce:)?source-text[^>]*>(.*?)</(?:ce:)?source-text>",
            content,
            re.DOTALL | re.IGNORECASE,
        )
    ]

    return {
        "title": title,
        "authors": deduplicate_preserve_order(authors),
        "keywords": deduplicate_preserve_order(keywords),
        "abstract": abstract,
        "body": body,
        "references": deduplicate_preserve_order(references),
    }


def parse_xml_root(input_file):
    """Parse XML and recover malformed documents when lxml is available."""
    try:
        return ET.parse(input_file).getroot(), "xml.etree"
    except ET.ParseError as exc:
        if LXML_ET is None:
            raise exc

        parser = LXML_ET.XMLParser(recover=True)
        root = LXML_ET.parse(str(input_file), parser).getroot()
        return root, "lxml-recover"


def extract_document_data(input_file):
    """Parse an XML file and dispatch to the best available extractor."""
    root, parser_name = parse_xml_root(input_file)
    tag = root.tag.lower()

    if "full-text-retrieval-response" in tag:
        print_status_block(
            "Detected Elsevier XML format based on root tag.",
            "Using Elsevier-specific extraction logic to handle embedded article content.",
        )
        return extract_elsevier_data(root), parser_name

    if tag.endswith("collection") or get_tag_name(root) == "collection":
        print_status_block(
            "Detected BioC XML format based on root tag.",
            "Using BioC-specific extraction logic to handle PubMed Central conversions and similar formats.",
        )
        return extract_bioc_data(root), parser_name

    article_root = find_elsevier_article_root(root)
    return {
        "title": extract_title(root, article_root),
        "authors": extract_authors(root, article_root),
        "keywords": extract_keywords(root, article_root),
        "abstract": extract_abstract(root, article_root),
        "body": extract_body(root, article_root),
        "references": extract_references(root, article_root),
    }, parser_name
def build_markdown(data, remove_references=False):
    """Render extracted publication content into a Markdown document."""
    parts = [f"# {data['title']}"]

    if data.get("authors"):
        parts.append(f"**Authors:** {', '.join(data['authors'])}")

    if data.get("keywords"):
        parts.append(f"**Keywords:** {', '.join(data['keywords'])}")

    if data.get("abstract"):
        parts.append("## Abstract")
        parts.append(data["abstract"])

    if data.get("body"):
        parts.append(data["body"])

    if not remove_references and data.get("references"):
        reference_lines = ["## References"]
        reference_lines.extend(f"{index}. {reference}" for index, reference in enumerate(data["references"], 1))
        parts.append("\n".join(reference_lines))

    return "\n\n".join(part for part in parts if part).strip() + "\n"


def resolve_output_path(input_path, output_path=None, output_dir=None):
    """Resolve where a Markdown file should be written for a given XML input."""
    input_path = Path(input_path)

    if output_path is not None:
        return Path(output_path)

    if output_dir is not None:
        return Path(output_dir) / f"{input_path.stem}.md"

    return input_path.with_suffix(".md")
def xml_to_md(filepath_xml, output_dir=None, remove_references=False, output_path=None):
    """Extract a single XML file into a Markdown file."""
    print_status_block("Parsing XML file.", f"Input: {filepath_xml}")

    try:
        data, parser_name = extract_document_data(filepath_xml)
        if parser_name == "lxml-recover":
            print_status_block(
                "Recovered malformed XML with lxml recovery parser.",
                f"Input: {filepath_xml}",
            )
    except ET.ParseError as exc:
        print_status_block(
            "Error parsing XML with xml.etree.",
            f"Input: {filepath_xml}",
            f"Details: {exc}",
            "Falling back to regex-based extraction due to malformed XML.",
        )
        data = extract_regex_fallback(filepath_xml)

    output_path = resolve_output_path(filepath_xml, output_path=output_path, output_dir=output_dir)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(build_markdown(data, remove_references=remove_references), encoding="utf-8")
    print_status_block("Successfully converted paper XML to Markdown.", f"Output: {output_path}")
    
    return output_path


def iter_input_files(input_path, recursive=False):
    """Yield XML input files from a file or directory input path."""
    input_path = Path(input_path)
    if input_path.is_file():
        return [input_path]
    if input_path.is_dir():
        pattern = "**/*.xml" if recursive else "*.xml"
        return sorted(input_path.glob(pattern))
    raise FileNotFoundError(f"Input path not found: {input_path}")


def main():
    parser = argparse.ArgumentParser(description="Convert XML article files into Markdown.")
    parser.add_argument(
        "input_path",
        nargs="?",
        default=None,
        help="Path to an XML file or a directory containing XML files.",
    )
    parser.add_argument(
        "-i",
        "--input-dir",
        help="Directory containing XML files to convert.",
    )
    parser.add_argument(
        "--output",
        help="Markdown output path for a single input XML file.",
    )
    parser.add_argument(
        "-d",
        "--output-dir",
        help="Directory for generated Markdown files. Recommended for directory inputs.",
    )
    parser.add_argument(
        "--recursive",
        action="store_true",
        help="Recursively search for XML files when the input path is a directory.",
    )
    parser.add_argument(
        "--no-references",
        action="store_true",
        help="Omit the references section from generated Markdown output.",
    )
    args = parser.parse_args()

    if args.input_dir and args.input_path:
        parser.error("Use either the positional input_path or --input-dir, not both.")

    input_source = args.input_dir or args.input_path or "10_1016_j_neuron_2025_05_018.xml"
    input_path = Path(input_source)
    if input_path.is_dir() and args.output:
        parser.error("--output can only be used with a single XML input file.")

    input_files = iter_input_files(input_path, recursive=args.recursive)
    if not input_files:
        print_status_block("No XML files found.", f"Search path: {input_path}")
        return

    for xml_file in input_files:
        xml_to_md(
            xml_file,
            output_dir=args.output_dir,
            remove_references=args.no_references,
            output_path=args.output if xml_file == input_files[0] and len(input_files) == 1 else None,
        )


if __name__ == "__main__":
    main()