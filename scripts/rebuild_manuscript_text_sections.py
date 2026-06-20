from pathlib import Path
import re

from docx import Document


ROOT = Path(__file__).resolve().parents[1]
DOCX = ROOT / "manuscript" / "Manuscript DNMT3L Neuron 20260617.docx"
SOURCES = [
    ROOT / "manuscript" / "Title_Abstract_DNMT3L_Neuron_20260617.md",
    ROOT / "manuscript" / "Introduction_section_DNMT3L_Neuron_20260617.md",
    ROOT / "manuscript" / "Results_section_DNMT3L_Neuron_20260615.md",
    ROOT / "manuscript" / "Discussion_section_DNMT3L_Neuron_20260617.md",
]
METHODS_MD = ROOT / "manuscript" / "Methods_section_DNMT3L_Neuron_20260617.md"


def split_markdown_blocks(text: str):
    blocks = []
    buf = []
    for line in text.splitlines():
        if not line.strip():
            if buf:
                blocks.append(" ".join(x.strip() for x in buf))
                buf = []
            continue
        buf.append(line)
    if buf:
        blocks.append(" ".join(x.strip() for x in buf))
    return blocks


def clean_literal_text(text: str) -> str:
    return text.replace("`", "")


def add_inline_runs(paragraph, text: str, *, bold=False, italic=False):
    text = clean_literal_text(text)
    pos = 0
    for match in re.finditer(r"\*([^*]+)\*|\^([0-9,*]+)", text):
        if match.start() > pos:
            run = paragraph.add_run(text[pos:match.start()])
            run.bold = bold
            run.italic = italic
        if match.group(1) is not None:
            run = paragraph.add_run(match.group(1))
            run.bold = bold
            run.italic = True
        else:
            run = paragraph.add_run(match.group(2))
            run.bold = bold
            run.italic = italic
            run.font.superscript = True
        pos = match.end()
    if pos < len(text):
        run = paragraph.add_run(text[pos:])
        run.bold = bold
        run.italic = italic


def insert_blocks_before(doc, anchor, blocks):
    for block in blocks:
        paragraph = doc.add_paragraph()
        if block.startswith("# "):
            add_inline_runs(paragraph, block[2:].strip(), bold=True)
        elif block == "## Abstract":
            add_inline_runs(paragraph, "Abstract", bold=True)
        elif block.startswith("## "):
            add_inline_runs(paragraph, block[3:].strip(), italic=True)
        elif block.startswith("*Correspondence:"):
            add_inline_runs(paragraph, block[1:].strip(), italic=True)
        else:
            add_inline_runs(paragraph, block)
        anchor._p.addprevious(paragraph._p)

    spacer = doc.add_paragraph()
    anchor._p.addprevious(spacer._p)


def find_paragraph_index(doc, text: str, start: int = 0) -> int:
    for i, paragraph in enumerate(doc.paragraphs[start:], start=start):
        if paragraph.text.strip() == text:
            return i
    raise RuntimeError(f"Could not find paragraph anchor: {text}")


def remove_paragraph_slice(doc, start: int, end: int):
    for paragraph in list(doc.paragraphs[start:end]):
        paragraph._element.getparent().remove(paragraph._element)


def main():
    doc = Document(DOCX)

    main_figures_idx = find_paragraph_index(doc, "Main Figures")
    main_figures_anchor = doc.paragraphs[main_figures_idx]
    remove_paragraph_slice(doc, 0, main_figures_idx)

    front_blocks = []
    for source in SOURCES:
        front_blocks.extend(split_markdown_blocks(source.read_text()))
    insert_blocks_before(doc, main_figures_anchor, front_blocks)

    methods_idx = find_paragraph_index(doc, "Methods")
    extended_idx = find_paragraph_index(doc, "Extended Data Figures", start=methods_idx + 1)
    extended_anchor = doc.paragraphs[extended_idx]
    remove_paragraph_slice(doc, methods_idx, extended_idx)
    methods_blocks = split_markdown_blocks(METHODS_MD.read_text())
    insert_blocks_before(doc, extended_anchor, methods_blocks)

    doc.save(DOCX)


if __name__ == "__main__":
    main()
