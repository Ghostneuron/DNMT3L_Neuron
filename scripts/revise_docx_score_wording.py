from pathlib import Path
import sys

from docx import Document


DEFAULT_DOCX = Path("manuscript/Manuscript DNMT3L Neuron 20260618_MolNeurobio_revised_v2.docx")


REPLACEMENTS = {
    "an integrated convergence score": "a heuristic integrated evidence-prioritization score",
    "integrated convergence score": "heuristic integrated evidence-prioritization score",
    "Top integrated convergence candidates": "Top-ranked evidence-prioritization candidates",
    "top integrated convergence candidates": "top-ranked evidence-prioritization candidates",
}

METHODS_ANCHOR = "Candidate scores were heuristic, capped prioritization scores rather than fitted model coefficients."
METHODS_PREFIX = (
    "The heuristic integrated evidence-prioritization score used for candidate ranking refers to "
    "the capped component-score framework below; it was not fit as a statistical model. "
)


def replace_in_runs(doc: Document) -> int:
    changed = 0
    for paragraph in doc.paragraphs:
        for run in paragraph.runs:
            new_text = run.text
            for old, new in REPLACEMENTS.items():
                new_text = new_text.replace(old, new)
            if new_text != run.text:
                run.text = new_text
                changed += 1
    return changed


def add_methods_bridge(doc: Document) -> bool:
    for paragraph in doc.paragraphs:
        if METHODS_ANCHOR in paragraph.text and METHODS_PREFIX.strip() not in paragraph.text:
            if paragraph.runs:
                paragraph.runs[0].text = METHODS_PREFIX + paragraph.runs[0].text
            else:
                paragraph.add_run(METHODS_PREFIX + paragraph.text)
            return True
    return False


def main() -> None:
    path = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_DOCX
    doc = Document(path)
    replacements = replace_in_runs(doc)
    methods_bridge_added = add_methods_bridge(doc)
    doc.save(path)
    print(f"Saved {path}")
    print(f"Phrase replacement runs touched: {replacements}")
    print(f"Methods bridge added: {methods_bridge_added}")


if __name__ == "__main__":
    main()
