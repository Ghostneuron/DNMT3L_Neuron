#!/usr/bin/env python3
"""Add fetal brain Hi-C TAD support wording to the formatted DOCX."""

from __future__ import annotations

from pathlib import Path

from docx import Document
from docx.oxml import OxmlElement
from docx.text.paragraph import Paragraph


ROOT = Path(__file__).resolve().parents[1]
DOCX = ROOT / "manuscript" / "Manuscript DNMT3L Neuron 20260618_MolNeurobio_revised_v2.docx"


RESULTS_LINK_PARAGRAPH = (
    "Because nearest-gene windows can over-assign methylation features to candidate genes, the final 24-gene set was re-evaluated with stricter link classes. "
    "All 24 final candidates retained promoter or gene-body support in at least one methylation layer, and none was supported only by distal 50-100 kb nearest-gene evidence. "
    "The five primary candidates retained strict human promoter or gene-body CpG support and mouse gene-body DMR support. "
    "Here, human triple-support CpGs were defined as human DNMT3L-overexpression-sensitive CpGs that also showed same-direction methylation in the SH-SY5Y phase-3 validation dataset and concordant methylation in at least one human DS methylation context. "
    "A fetal-brain Hi-C domain-level sensitivity analysis further showed that 19 of 20 candidates with human DNMT3L-sensitive CpG-gene links had at least one CpG and candidate-gene TSS within the same Won et al. fetal cortical plate or germinal-zone TAD [39]; 16 candidates retained same-TAD support among strict fetal-positive human probes. "
    "This result supports fetal cortical regulatory-domain plausibility but does not establish direct enhancer-promoter contact. "
    "Across the broader set, 17 of 24 candidates had human triple-support CpGs and all 24 had mouse Dnmt3l-sensitive DMR support; however, the expression/pathway tier lacked human triple-support CpGs and should therefore be read as mouse/context-supported and methylation-unresolved rather than as a human methylation-convergent tier (Supplementary Tables 6d, 8h, and 8i)."
)


RESULTS_REG_PARAGRAPH = (
    "The primary candidates were not selected because they encode methyltransferases or STAT proteins. "
    "Instead, they emerged as developmental regulators or signaling genes near DNMT3L-linked methylation features with cross-context regulatory support. "
    "Candidate-level support was retained as separate evidence components, including human methylation, human expression, mouse evidence, cross-species expression priority, and regulatory-context scores, so the final rank should be read as prioritization rather than as a causal weight. "
    "Candidate-linked human DNMT3L-sensitive CpG-gene records were enriched for fetal brain regulatory annotations after matched permutation testing. "
    "In the primary-five set, human DNMT3L-linked CpG-gene records overlapped ATAC/open-chromatin annotations above expectation (12 of 41 records; 2.44-fold; empirical P = 0.0025), H3K27ac (27 of 41; 2.28-fold; empirical P = 0.0001), and H3K27me3 (28 of 41; 2.71-fold; empirical P = 0.0001). "
    "The broader 24-gene set also showed regulatory enrichment, including strict human triple-support probes (54 of 55; 1.39-fold; empirical P = 0.0001). "
    "The fetal Hi-C TAD analysis placed all five primary candidates in same-domain CP/GZ contexts for their strict human support probes. "
    "Mouse Dnmt3l-sensitive DMRs supported individual loci but were not globally enriched above matched mouse DMR background for ATAC or histone-mark overlap, indicating that the mouse evidence was stronger as locus-level support than as a broad annotation-enrichment result (Fig. 2; Extended Data Fig. 5; Supplementary Tables 6 and 8)."
)


METHODS_EVIDENCE_PARAGRAPH = (
    "The integrated candidate evidence was treated as a prioritization framework rather than as a fitted causal model. "
    "Component scores were retained separately in Supplementary Table 6b, including mouse evidence score, human methylation support score, human expression support score, cross-species expression priority score, and regulatory-context score. "
    "Final tier labels were assigned from combined evidence class and manual review of the final candidate set, not from a single unvalidated mathematical weight intended to imply mechanism. "
    "To add a fetal-brain three-dimensional-genome sensitivity layer, GSE157078 hg19 CpG coordinates and GRCh37 candidate-gene transcription start sites were intersected with Won et al. GSE77565 fetal cortical plate and germinal-zone Hi-C TADs [39]. "
    "CpG-gene links were annotated as same-TAD supported when the CpG and candidate-gene TSS fell within the same CP or GZ TAD. "
    "This analysis was used as domain-level regulatory plausibility, not as evidence for direct enhancer-promoter looping; direct chromatin contacts, ABC-style enhancer-gene scores, or perturbation data remain necessary for causal validation."
)


DISCUSSION_SUPPORT_PARAGRAPH = (
    "Several lines of evidence support a methylation-linked mechanism. "
    "First, DNMT3L/Dnmt3l overexpression produced measurable, mostly hypermethylation-skewed methylation changes in human fetal neural and mouse postnatal cortical contexts. "
    "Second, DS brain methylomes showed partial concordance with DNMT3L-sensitive methylation, strongest in selected fetal, adult bulk/glial, and Dp10 contexts. "
    "Third, integrated methylation, expression, and regulatory analyses nominated 24 downstream candidates, including primary candidates PLPP3, FOXP1, FGFR2, ZFHX3, and SIAH3, and 16 of these 24 genes have neurodevelopmental relevance. "
    "This candidate set was not dependent on distal 100 kb-only assignment: all 24 retained promoter/gene-body support in at least one methylation layer, and all five primary candidates retained strict human promoter/gene-body CpG support plus mouse gene-body DMR support. "
    "Fourth, candidate-linked human DNMT3L-sensitive CpG-gene records were enriched for fetal brain regulatory annotations, including open chromatin and H3K27ac/H3K27me3, and most human CpG-gene candidate links fell within the same fetal cortical plate or germinal-zone Hi-C TAD as the candidate-gene TSS. "
    "Fifth, in the mouse Dnmt3l-OE postnatal cortex, glial/myelination and neurodevelopmental signaling modules were enriched for nearby Dnmt3l-sensitive DMRs and for DNMT3A-linked methylation-context overlap. "
    "Together, these findings indicate that DNMT3L dosage can be associated with altered methylation at developmental regulatory loci."
)


DISCUSSION_LIMIT_PARAGRAPH = (
    "The evidence also argues against overextending the methylation model. "
    "DNMT3L-sensitive methylation did not fully reproduce DS methylation, and DS methylation patterns were region-, cell-type-, and stage-dependent. "
    "Most candidate gene links are based on proximity, regulatory annotation, same-TAD domain support, or cross-context concordance, not on paired methylation-expression causality in the same cells. "
    "Same fetal Hi-C TAD support strengthens regulatory-domain plausibility but is weaker than direct enhancer-promoter contact. "
    "This is especially important for directionally complex loci. "
    "FOXP1 and parts of the SIAH3 evidence fit a repression-like interpretation in at least one context, whereas PLPP3, FGFR2, and ZFHX3 more often show methylation-expression relationships compatible with enhancer methylation, gene-body methylation, indirect timing effects, or cell-state shifts rather than simple promoter silencing. "
    "The available Dnmt3a-cKO transcriptome did not recapitulate Dnmt3l-OE expression effects, but that does not prove DNMT3A independence because the datasets differ in age, cell type, perturbation direction, and assay context. "
    "Likewise, HOMER and motif analyses identify sequence-context hypotheses rather than direct transcription-factor binding or methylation-mediated regulation. "
    "The strongest statement is therefore that methylation-linked regulation is the better-supported arm of the current in-silico evidence, not that every downstream expression change is directly caused by local DNA methylation."
)


FUTURE_DIRECTIONS_PARAGRAPH = (
    "Second, direct DNMT3A/DNMT3B dependency should be tested in matched neural systems. "
    "DNMT3L perturbation combined with DNMT3A or DNMT3B knockdown, catalytic inhibition, or rescue would clarify whether candidate methylation changes require de novo methyltransferase activity. "
    "CUT&RUN, CUT&Tag, or ChIP-seq for DNMT3A, DNMT3B, H3K27me3, H3K27ac, and relevant developmental TFs at candidate loci would help move beyond proximity-based links. "
    "The fetal Hi-C TAD analysis added here helps identify candidate links that are at least compatible with fetal cortical three-dimensional domains, but promoter-capture Hi-C, high-resolution loop calls, or ABC-style enhancer-gene maps would still be needed to test whether candidate CpGs or DMRs physically contact candidate promoters. "
    "Targeted bisulfite or long-read methylation assays across candidate regulatory regions could validate locus-level methylation architecture."
)


WON_REFERENCE = (
    "39. Won H, de la Torre-Ubieta L, Stein JL, Parikshak NN, Huang J, Opland CK, Gandal MJ, Sutton GJ, Hormozdiari F, Lu D, Lee C, Eskin E, Voineagu I, Ernst J, Geschwind DH (2016) "
    "Chromosome conformation elucidates regulatory relationships in developing human brain. Nature 538 (7626):523-527. doi:10.1038/nature19847"
)


def replace_paragraph_text(paragraph, text: str) -> None:
    paragraph.clear()
    paragraph.add_run(text)


def insert_after(paragraph, text: str) -> Paragraph:
    new_p = OxmlElement("w:p")
    paragraph._p.addnext(new_p)
    new_para = Paragraph(new_p, paragraph._parent)
    new_para.style = paragraph.style
    new_para.add_run(text)
    return new_para


def main() -> None:
    doc = Document(DOCX)
    hits = {
        "dataset_list": False,
        "results_links": False,
        "results_reg": False,
        "methods": False,
        "discussion_support": False,
        "discussion_limit": False,
        "future": False,
        "reference": False,
    }

    for para in doc.paragraphs:
        text = para.text
        if "The main dataset used were:" in text and "GSE149268 [28]; GSE162170" in text:
            replace_paragraph_text(
                para,
                text.replace("GSE149268 [28]; GSE162170", "GSE149268 [28]; GSE77565 fetal brain Hi-C [39]; GSE162170"),
            )
            hits["dataset_list"] = True
        elif text.startswith("Because nearest-gene windows can over-assign methylation features"):
            replace_paragraph_text(para, RESULTS_LINK_PARAGRAPH)
            hits["results_links"] = True
        elif text.startswith("The primary candidates were not selected because they encode methyltransferases"):
            replace_paragraph_text(para, RESULTS_REG_PARAGRAPH)
            hits["results_reg"] = True
        elif text.startswith("The integrated candidate evidence was treated as a prioritization framework"):
            replace_paragraph_text(para, METHODS_EVIDENCE_PARAGRAPH)
            hits["methods"] = True
        elif text.startswith("Several lines of evidence support a methylation-linked mechanism"):
            replace_paragraph_text(para, DISCUSSION_SUPPORT_PARAGRAPH)
            hits["discussion_support"] = True
        elif text.startswith("The evidence also argues against overextending the methylation model"):
            replace_paragraph_text(para, DISCUSSION_LIMIT_PARAGRAPH)
            hits["discussion_limit"] = True
        elif text.startswith("Second, direct DNMT3A/DNMT3B dependency should be tested"):
            replace_paragraph_text(para, FUTURE_DIRECTIONS_PARAGRAPH)
            hits["future"] = True
        elif text.startswith("39. Won H"):
            replace_paragraph_text(para, WON_REFERENCE)
            hits["reference"] = True

    if not hits["reference"]:
        for para in doc.paragraphs:
            if para.text.startswith("38. Heyn P"):
                insert_after(para, WON_REFERENCE)
                hits["reference"] = True
                break

    missing = [key for key, ok in hits.items() if not ok]
    if missing:
        raise RuntimeError(f"Missing DOCX edit anchors: {missing}")

    doc.save(DOCX)
    print(f"Saved {DOCX}")
    print("Updated fetal Hi-C/TAD support wording and added Won et al. reference 39.")


if __name__ == "__main__":
    main()
