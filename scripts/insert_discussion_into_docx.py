from pathlib import Path
import re

from docx import Document


ROOT = Path(__file__).resolve().parents[1]
DOCX = ROOT / "manuscript" / "Manuscript DNMT3L Neuron 20260617.docx"
DISCUSSION_MD = ROOT / "manuscript" / "Discussion_section_DNMT3L_Neuron_20260617.md"


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


def add_markdown_runs(paragraph, text: str):
    pos = 0
    for match in re.finditer(r"\*([^*]+)\*", text):
        if match.start() > pos:
            paragraph.add_run(text[pos:match.start()])
        run = paragraph.add_run(match.group(1))
        run.italic = True
        pos = match.end()
    if pos < len(text):
        paragraph.add_run(text[pos:])


def remove_existing_discussion(doc):
    discussion_idx = None
    main_figures_idx = None
    for i, paragraph in enumerate(doc.paragraphs):
        text = paragraph.text.strip()
        if text == "Discussion" and discussion_idx is None:
            discussion_idx = i
        if text == "Main Figures" and main_figures_idx is None:
            main_figures_idx = i
    if main_figures_idx is None:
        raise RuntimeError("Could not find Main Figures section.")
    if discussion_idx is not None and discussion_idx < main_figures_idx:
        for paragraph in list(doc.paragraphs[discussion_idx:main_figures_idx]):
            paragraph._element.getparent().remove(paragraph._element)


def find_anchor(doc, anchor_text):
    for paragraph in doc.paragraphs:
        if paragraph.text.strip() == anchor_text:
            return paragraph
    raise RuntimeError(f"Could not find anchor: {anchor_text}")


def main():
    doc = Document(DOCX)
    remove_existing_discussion(doc)
    anchor = find_anchor(doc, "Main Figures")

    for block in split_markdown_paragraphs(DISCUSSION_MD.read_text()):
        paragraph = doc.add_paragraph()
        if block.startswith("# "):
            run = paragraph.add_run(block[2:].strip())
            run.bold = True
        elif block.startswith("## "):
            run = paragraph.add_run(block[3:].strip())
            run.italic = True
        else:
            add_markdown_runs(paragraph, block)
        anchor._p.addprevious(paragraph._p)

    spacer = doc.add_paragraph()
    anchor._p.addprevious(spacer._p)
    doc.save(DOCX)


if __name__ == "__main__":
    main()
