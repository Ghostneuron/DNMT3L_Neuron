from pathlib import Path
import re

from docx import Document


ROOT = Path(__file__).resolve().parents[1]
DOCX = ROOT / "manuscript" / "Manuscript DNMT3L Neuron 20260617.docx"
TITLE_ABSTRACT_MD = ROOT / "manuscript" / "Title_Abstract_DNMT3L_Neuron_20260617.md"


def split_markdown_paragraphs(text: str):
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


def add_markdown_runs(paragraph, text: str, *, bold=False):
    pos = 0
    for match in re.finditer(r"\*([^*]+)\*", text):
        if match.start() > pos:
            run = paragraph.add_run(text[pos:match.start()])
            run.bold = bold
        run = paragraph.add_run(match.group(1))
        run.bold = bold
        run.italic = True
        pos = match.end()
    if pos < len(text):
        run = paragraph.add_run(text[pos:])
        run.bold = bold


def find_anchor(doc, anchor_text):
    for paragraph in doc.paragraphs:
        if paragraph.text.strip() == anchor_text:
            return paragraph
    raise RuntimeError(f"Could not find anchor: {anchor_text}")


def remove_existing_front_matter(doc):
    intro_idx = None
    for i, paragraph in enumerate(doc.paragraphs):
        if paragraph.text.strip() == "Introduction":
            intro_idx = i
            break
    if intro_idx is None:
        raise RuntimeError("Could not find Introduction section.")
    for paragraph in list(doc.paragraphs[:intro_idx]):
        paragraph._element.getparent().remove(paragraph._element)


def main():
    doc = Document(DOCX)
    remove_existing_front_matter(doc)
    anchor = find_anchor(doc, "Introduction")

    for block in split_markdown_paragraphs(TITLE_ABSTRACT_MD.read_text()):
        paragraph = doc.add_paragraph()
        if block.startswith("# "):
            add_markdown_runs(paragraph, block[2:].strip(), bold=True)
        elif block.startswith("## "):
            run = paragraph.add_run(block[3:].strip())
            run.bold = True
        else:
            add_markdown_runs(paragraph, block)
        anchor._p.addprevious(paragraph._p)

    spacer = doc.add_paragraph()
    anchor._p.addprevious(spacer._p)
    doc.save(DOCX)


if __name__ == "__main__":
    main()
