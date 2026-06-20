#!/usr/bin/env python3
"""Submission-cleanup pass for the formatted Molecular Neurobiology DOCX.

This script edits a working DOCX copy in place while preserving the EndNote
bibliography paragraphs. It removes confusing "Corrected" wording, expands
STAT1/3 shorthand to explicit STAT1/STAT3, and italicizes curated gene symbols
in manuscript text and figure legends before the bibliography.
"""

from __future__ import annotations

from copy import deepcopy
from pathlib import Path
import re
import sys

from docx import Document
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.text.paragraph import Paragraph


ROOT = Path(__file__).resolve().parents[1]
DOCX = (
    Path(sys.argv[1]).resolve()
    if len(sys.argv) > 1
    else ROOT / "manuscript" / "Manuscript DNMT3L Neuron 20260618_MolNeurobio_revised.docx"
)

GENE_SYMBOLS = [
    "DNMT3L",
    "DNMT3A",
    "DNMT3B",
    "Dnmt3l",
    "Dnmt3a",
    "Dnmt3b",
    "STAT1",
    "STAT3",
    "Stat1",
    "Stat3",
    "HES3",
    "ASCL1",
    "NEUROD2",
    "NEUROG2",
    "CDKN1A",
    "PLPP3",
    "FOXP1",
    "FGFR2",
    "ZFHX3",
    "SIAH3",
    "FAM222A",
    "COL26A1",
    "KIRREL3",
    "CAMTA1",
    "LMX1B",
    "OTX1",
    "OTX2",
    "LHX2",
    "DRAXIN",
    "NFIA",
    "PRDM16",
    "NR2F1",
    "FZD7",
    "TGIF2",
    "PRDM8",
    "CRABP2",
    "BCL11B",
    "RXRA",
    "CYP26A1",
    "Plpp3",
    "Foxp1",
    "Fgfr2",
    "Zfhx3",
    "Siah3",
    "Fam222a",
    "Col26a1",
    "Kirrel3",
    "Camta1",
    "Lmx1b",
    "Otx1",
    "Otx2",
    "Lhx2",
    "Draxin",
    "Nfia",
    "Prdm16",
    "Nr2f1",
    "Fzd7",
    "Tgif2",
    "Prdm8",
    "Crabp2",
    "Bcl11b",
    "Rxra",
    "Cyp26a1",
    "CUX2",
    "TLE4",
    "SST",
    "RORB",
    "RELN",
    "CALB2",
]

GENE_PATTERN = re.compile(
    r"(?<![A-Za-z0-9])(" + "|".join(re.escape(s) for s in sorted(GENE_SYMBOLS, key=len, reverse=True)) + r")(?![A-Za-z0-9])"
)

TEXT_REPLACEMENTS = {
    "Corrected fetal bulk DS expression analyses": "Fetal bulk DS expression analyses",
    "Corrected human DNMT3L-OE count-based expression values": "Human DNMT3L-OE count-based expression values",
    "DNMT3L-STAT1/3": "DNMT3L-STAT1/STAT3",
    "STAT1/3 phosphorylation": "STAT1/STAT3 phosphorylation",
    "DNMT3L-STAT1/3 differentiation axis": "DNMT3L-STAT1/STAT3 differentiation axis",
    "STAT1/3 differentiation axis": "STAT1/STAT3 differentiation axis",
    "STAT1/3 activity": "STAT1/STAT3 activity",
}


def paragraph_after(paragraph, text: str | None = None) -> Paragraph:
    new_p = OxmlElement("w:p")
    paragraph._p.addnext(new_p)
    new_para = Paragraph(new_p, paragraph._parent)
    if text:
        new_para.add_run(text)
    return new_para


def paragraph_before(paragraph, text: str | None = None) -> Paragraph:
    new_p = OxmlElement("w:p")
    paragraph._p.addprevious(new_p)
    new_para = Paragraph(new_p, paragraph._parent)
    if text:
        new_para.add_run(text)
    return new_para


def add_bold_heading(paragraph: Paragraph, text: str) -> None:
    run = paragraph.add_run(text)
    run.bold = True


def add_keywords_if_missing(doc: Document) -> None:
    if any(p.text.strip().startswith("Keywords:") for p in doc.paragraphs):
        return
    for i, para in enumerate(doc.paragraphs):
        if para.text.strip() == "Abstract" and i + 1 < len(doc.paragraphs):
            kw = paragraph_after(doc.paragraphs[i + 1])
            run = kw.add_run("Keywords: ")
            run.bold = True
            kw.add_run("Down syndrome; ")
            gene = kw.add_run("DNMT3L")
            gene.italic = True
            kw.add_run("; DNA methylation; neurodevelopment; public multi-omics; STAT signaling")
            return


def add_declaration_placeholders_if_missing(doc: Document) -> None:
    if any(p.text.strip() == "Statements and Declarations" for p in doc.paragraphs):
        return

    data_para = next((p for p in doc.paragraphs if p.text.strip() == "Data availability"), None)
    if data_para is None:
        return

    heading = paragraph_before(data_para)
    add_bold_heading(heading, "Statements and Declarations")

    funding_h = paragraph_before(data_para)
    add_bold_heading(funding_h, "Funding")
    paragraph_before(data_para, "No external funding was received for this study.")

    ethics_h = paragraph_before(data_para)
    add_bold_heading(ethics_h, "Ethics approval and consent to participate")
    paragraph_before(
        data_para,
        "This study analyzed public datasets, public source tables, or public visualization outputs and generated no new human or animal data. Additional ethics approval and participant consent were not required for this in-silico analysis.",
    )

    consent_h = paragraph_before(data_para)
    add_bold_heading(consent_h, "Consent for publication")
    paragraph_before(data_para, "Not applicable.")


def clone_run(base_r, text: str, italic: bool) -> OxmlElement:
    new_r = OxmlElement("w:r")
    rpr = base_r.rPr
    if rpr is not None:
        new_r.append(deepcopy(rpr))

    if italic:
        new_r_rpr = new_r.find(qn("w:rPr"))
        if new_r_rpr is None:
            new_r_rpr = OxmlElement("w:rPr")
            new_r.insert(0, new_r_rpr)
        for tag in ("w:i", "w:iCs"):
            for existing in list(new_r_rpr.findall(qn(tag))):
                new_r_rpr.remove(existing)
            el = OxmlElement(tag)
            el.set(qn("w:val"), "1")
            new_r_rpr.append(el)

    t = OxmlElement("w:t")
    if text[:1].isspace() or text[-1:].isspace():
        t.set(qn("xml:space"), "preserve")
    t.text = text
    new_r.append(t)
    return new_r


def apply_run_replacements(paragraph) -> int:
    n = 0
    for run in paragraph.runs:
        text = run.text
        new_text = text
        for old, new in TEXT_REPLACEMENTS.items():
            new_text = new_text.replace(old, new)
        if new_text != text:
            run.text = new_text
            n += 1
    return n


def italicize_gene_runs(paragraph) -> int:
    n = 0
    for run in list(paragraph.runs):
        text = run.text
        if not text or not GENE_PATTERN.search(text):
            continue

        pieces: list[tuple[str, bool]] = []
        last = 0
        for match in GENE_PATTERN.finditer(text):
            if match.start() > last:
                pieces.append((text[last : match.start()], False))
            pieces.append((match.group(1), True))
            n += 1
            last = match.end()
        if last < len(text):
            pieces.append((text[last:], False))

        parent = run._r.getparent()
        idx = parent.index(run._r)
        for offset, (piece, is_gene) in enumerate(pieces):
            if piece:
                parent.insert(idx + offset, clone_run(run._r, piece, is_gene))
        parent.remove(run._r)
    return n


def iter_editable_paragraphs(doc: Document):
    for para in doc.paragraphs:
        if para.style.name == "EndNote Bibliography":
            break
        yield para


def main() -> None:
    doc = Document(DOCX)
    add_keywords_if_missing(doc)
    add_declaration_placeholders_if_missing(doc)
    replacement_count = 0
    italic_count = 0
    for para in iter_editable_paragraphs(doc):
        replacement_count += apply_run_replacements(para)
        italic_count += italicize_gene_runs(para)
    doc.save(DOCX)
    print(f"Saved {DOCX}")
    print(f"Text replacement runs touched: {replacement_count}")
    print(f"Gene-symbol italic matches: {italic_count}")


if __name__ == "__main__":
    main()
