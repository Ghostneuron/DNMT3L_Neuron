#!/usr/bin/env python3
"""Address final manuscript follow-up points in the formatted DOCX copy."""

from __future__ import annotations

from pathlib import Path

from docx import Document
from docx.oxml import OxmlElement
from docx.text.paragraph import Paragraph


ROOT = Path(__file__).resolve().parents[1]
DOCX = ROOT / "manuscript" / "Manuscript DNMT3L Neuron 20260618_MolNeurobio_revised_v2.docx"


TRIPLE_SUPPORT_PARAGRAPH = (
    "Because nearest-gene windows can over-assign methylation features to candidate genes, the final 24-gene set was re-evaluated with stricter link classes. "
    "All 24 final candidates retained promoter or gene-body support in at least one methylation layer, and none was supported only by distal 50-100 kb nearest-gene evidence. "
    "The five primary candidates retained strict human promoter or gene-body CpG support and mouse gene-body DMR support. "
    "Here, human triple-support CpGs were defined as human DNMT3L-overexpression-sensitive CpGs that also showed same-direction methylation in the SH-SY5Y phase-3 validation dataset and concordant methylation in at least one human DS methylation context. "
    "Across the broader set, 17 of 24 candidates had human triple-support CpGs and all 24 had mouse Dnmt3l-sensitive DMR support; however, the expression/pathway tier lacked human triple-support CpGs and should therefore be read as mouse/context-supported and methylation-unresolved rather than as a human methylation-convergent tier (Supplementary Tables 8h and 8i)."
)


SCORING_METHODS_PARAGRAPH = (
    "The heuristic integrated evidence-prioritization score used for candidate ranking refers to the capped component-score framework below; it was not fit as a statistical model. "
    "Candidate scores were heuristic, capped prioritization scores rather than fitted model coefficients. "
    "The mouse evidence score summed a methylation relationship score, an expression score, and a Dp10 support score. "
    "The methylation relationship score counted promoter DMRs as 3 points each capped at 9, gene-body DMRs as 2 points each capped at 4, and nearest-100 kb DMRs as 1 point each capped at 2. "
    "The mouse expression score awarded one point each for detectable expression, absolute log2FC >= 0.58, absolute log2FC >= 1, FDR <= 0.05 with absolute log2FC >= 1, and sex-consistent direction; the Dp10 score awarded up to two points for detectable direction/effect and concordance with mouse Dnmt3l-OE direction. "
    "The human methylation support score was calculated as min(10, 3 x fetal mean-positive triple-support CpGs) + min(5, 2 x fetal-majority triple-support CpGs) + min(3, all human triple-support CpGs) + min(3, promoter triple-support CpGs) + min(2, gene-body triple-support CpGs). "
    "The human expression support score was min(2, same-direction bulk contexts) + min(4, same-direction fetal cell-type contexts), with one additional point when both bulk and cell-type support were present. "
    "The cross-species expression priority score equaled mouse evidence score plus human expression support score, and the combined methylation-expression score equaled human methylation support score plus cross-species expression priority score. "
    "Regulatory-context score added 3 points for strict human regulatory overlap or 1 point for broader human regulatory overlap, 2 points for mouse regulatory DMR overlap, 1 point when at least half of linked mouse DMRs overlapped mouse ATAC, 1 point for mouse H3K27ac overlap, and 1 point for mouse H3K27me3 overlap. "
    "These component scores were retained separately in Supplementary Table 6b; final tier labels were assigned from evidence class and manual review, not from a single causal score."
)


PRIMARY_CANDIDATE_DISCUSSION = (
    "The five primary candidates also represent biologically distinct neurodevelopmental entry points. "
    "FOXP1 is the clearest prior neurodevelopmental-disorder gene in this set, with relevance to cortical and striatal circuit development, speech/language phenotypes, intellectual disability, and autism-spectrum presentations; its repression-like methylation-expression pattern therefore has direct neurological plausibility. "
    "FGFR2 connects the candidate set to FGF-dependent progenitor proliferation, cortical patterning, and developmental syndrome biology. "
    "PLPP3 encodes a lipid-phosphate phosphatase and is more naturally interpreted through extracellular lipid signaling, cell migration, vascular-glial interfaces, and regional or cell-state regulation than through a simple neuronal transcription-factor model. "
    "ZFHX3 encodes a large homeobox/zinc-finger transcription factor connected to developmental transcriptional regulation and neuronal or circadian programs; in this study it is valuable mainly as a regulatory-convergent and relatively novel candidate because its mouse expression effect is weak after public count reanalysis. "
    "SIAH3 is the least established neurodevelopmental-disorder gene among the primary five, but its E3-ligase-family context suggests possible ubiquitin-linked, proteostasis, or differentiation-state mechanisms. "
    "Thus, the primary set should not be framed as a single downstream pathway; it spans transcriptional identity, growth-factor signaling, lipid/vascular-glial signaling, transcriptional timing, and ubiquitin-linked regulation."
)


def paragraph_after(paragraph, text: str) -> Paragraph:
    new_p = OxmlElement("w:p")
    paragraph._p.addnext(new_p)
    new_para = Paragraph(new_p, paragraph._parent)
    new_para.add_run(text)
    return new_para


def replace_paragraph_text(paragraph, text: str) -> None:
    paragraph.clear()
    paragraph.add_run(text)


def main() -> None:
    doc = Document(DOCX)

    replaced_triple = False
    inserted_scoring = False
    inserted_primary_discussion = False

    for i, para in enumerate(list(doc.paragraphs)):
        txt = para.text
        if "Across the broader set, 17 of 24 candidates had strict human triple-support CpGs" in txt:
            replace_paragraph_text(para, TRIPLE_SUPPORT_PARAGRAPH)
            replaced_triple = True

        if (
            "The integrated candidate evidence was treated as a prioritization framework"
            in txt
            and "Component scores were retained separately" in txt
            and "Supplementary Table 6b" in txt
        ):
            paragraph_after(para, SCORING_METHODS_PARAGRAPH)
            inserted_scoring = True

        if (
            "The evidence also argues against overextending the methylation model."
            in txt
            and "The strongest statement is therefore that methylation-linked regulation"
            in txt
        ):
            paragraph_after(para, PRIMARY_CANDIDATE_DISCUSSION)
            inserted_primary_discussion = True

    if not replaced_triple:
        raise RuntimeError("Could not find triple-support Results paragraph.")
    if not inserted_scoring:
        raise RuntimeError("Could not find Methods scoring insertion point.")
    if not inserted_primary_discussion:
        raise RuntimeError("Could not find Discussion candidate insertion point.")

    doc.save(DOCX)
    print(f"Saved {DOCX}")


if __name__ == "__main__":
    main()
