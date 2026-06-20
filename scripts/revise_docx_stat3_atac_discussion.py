from pathlib import Path
import sys

from docx import Document
from docx.oxml import OxmlElement
from docx.text.paragraph import Paragraph


DEFAULT_DOCX = Path("manuscript/Manuscript DNMT3L Neuron 20260618_MolNeurobio_revised_v2.docx")

ANCHOR = "However, public DS chromatin data did not support STAT as a dominant direct motif mechanism."

NEW_PARAGRAPH = (
    "The absence of a local ATAC signal at STAT3 should therefore not be interpreted as absence of STAT3 regulation. "
    "ATAC detects changes in chromatin accessibility, whereas STAT3 RNA or activity can increase from an already "
    "permissive promoter/enhancer state, distal regulatory elements outside the tested cis-window, cytokine/JAK-STAT "
    "feedback, phosphorylation-dependent pathway activation, or shifts toward glial, immune-reactive, or stress-response "
    "cell states. This distinction is important because DS expression and DNMT3L-overexpression neuronal contexts support "
    "recurrent STAT3 upregulation, whereas the chromatin analyses do not identify STAT3 as a local chromatin-opening target. "
    "In the present evidence map, STAT3 is therefore best interpreted as a signaling or cell-state effector that can be "
    "activated downstream of upstream pathway perturbation, rather than as a candidate whose upregulation requires detectable "
    "local ATAC opening."
)


def insert_after(paragraph, text: str) -> Paragraph:
    new_p = OxmlElement("w:p")
    paragraph._p.addnext(new_p)
    new_paragraph = Paragraph(new_p, paragraph._parent)
    new_paragraph.style = paragraph.style
    new_paragraph.add_run(text)
    return new_paragraph


def main() -> None:
    path = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_DOCX
    doc = Document(path)
    full_text = "\n".join(par.text for par in doc.paragraphs)
    if NEW_PARAGRAPH[:90] in full_text:
        print(f"Already present in {path}")
        return
    for paragraph in doc.paragraphs:
        if ANCHOR in paragraph.text:
            insert_after(paragraph, NEW_PARAGRAPH)
            doc.save(path)
            print(f"Saved {path}")
            print("Inserted STAT3/ATAC discussion paragraph")
            return
    raise SystemExit(f"Anchor paragraph not found in {path}")


if __name__ == "__main__":
    main()
