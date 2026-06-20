from pathlib import Path
import re

from docx import Document
from docx.enum.text import WD_BREAK


ROOT = Path(__file__).resolve().parents[1]
DOCX = ROOT / "manuscript" / "Manuscript DNMT3L Neuron 20260617.docx"
INTRO_MD = ROOT / "manuscript" / "Introduction_section_DNMT3L_Neuron_20260617.md"


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
    """Add runs, treating *...* spans as italics."""
    pos = 0
    for match in re.finditer(r"\*([^*]+)\*", text):
        if match.start() > pos:
            paragraph.add_run(text[pos:match.start()])
        run = paragraph.add_run(match.group(1))
        run.italic = True
        pos = match.end()
    if pos < len(text):
        paragraph.add_run(text[pos:])


def set_body_format(paragraph):
    paragraph.paragraph_format.space_after = None
    paragraph.paragraph_format.space_before = None
    paragraph.paragraph_format.line_spacing = None


def main():
    doc = Document(DOCX)

    # Remove an older inserted Introduction if this script is rerun.
    first_results_idx = None
    for i, paragraph in enumerate(doc.paragraphs):
        if paragraph.text.strip() == "Results":
            first_results_idx = i
            break
    if first_results_idx is None:
        raise RuntimeError("Could not find the Results section.")

    if first_results_idx > 0 and doc.paragraphs[0].text.strip() == "Introduction":
        for paragraph in list(doc.paragraphs[:first_results_idx]):
            paragraph._element.getparent().remove(paragraph._element)

    first_results = doc.paragraphs[0]
    blocks = split_markdown_paragraphs(INTRO_MD.read_text())

    for block in blocks:
        new_para = doc.add_paragraph()
        if block.startswith("# "):
            run = new_para.add_run(block[2:].strip())
            run.bold = True
        else:
            add_markdown_runs(new_para, block)
            set_body_format(new_para)
        first_results._p.addprevious(new_para._p)

    spacer = doc.add_paragraph()
    first_results._p.addprevious(spacer._p)

    doc.save(DOCX)


if __name__ == "__main__":
    main()
