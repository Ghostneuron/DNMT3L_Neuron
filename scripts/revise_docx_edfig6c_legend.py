#!/usr/bin/env python3
"""Revise active DOCX to treat fetal Hi-C TAD support as ED Fig. 6c."""

from __future__ import annotations

from pathlib import Path

from docx import Document


ROOT = Path(__file__).resolve().parents[1]
DOCX = ROOT / "manuscript" / "Manuscript DNMT3L Neuron 20260618_MolNeurobio_revised_v2.docx"


ED6_TITLE = "Extended Data Figure 6. Neurodevelopmental candidate annotation, Dnmt3a threshold sensitivity, and fetal Hi-C TAD support."


ED6_LEGEND = (
    "Functional annotation, Dnmt3a-cKO threshold-sensitivity, and fetal Hi-C domain-support analyses test whether final candidates are neurodevelopmentally relevant, mirrored by the available Dnmt3a loss-of-function transcriptome, and compatible with fetal cortical 3D regulatory domains. "
    "(a) Final-candidate annotation heatmap maps the 24 genes to curated neurodevelopmental categories: core neurogenesis, brain-region patterning, neurite/synapse/migration, glial context, and neurogenesis signaling. Row markers denote candidate tier: red, primary methylation-convergent; yellow, secondary methylation/expression; purple, developmental/Polycomb; green, expression/pathway methylation-unresolved. The Dnmt3l OE DE column uses the public GSE268951 expected-count DESeq2/edgeR reanalysis under the manuscript rule (FDR <= 0.05 and absolute log2FC >= 1). The Dnmt3a cKO DE column uses the public GSE141587 Dnmt3a-cKO edgeR supplementary table with the same FDR/logFC cutoff. Sixteen of 24 candidates are annotated as neurogenesis/neurodevelopment-related, whereas none of the final candidates passes the Dnmt3a-cKO DE threshold. "
    "(b) Dnmt3a-cKO threshold-sensitivity plot counts public GSE141587 Dnmt3a-cKO genes passing FDR <= 0.05 as the absolute logFC cutoff is relaxed from 1.0 to 0.58, 0.25, and 0.0. Lines show protein-coding Dnmt3a-cKO DE genes, neurogenesis-related Dnmt3a-cKO DE genes, and overlap with a Dnmt3l-OE annotation flag. The overlap line is not the main manuscript Dnmt3l-OE count-based DE result and is not derived from the unpublished worktable; for this robustness visualization only, it uses the local GSE268951 nominal flag in Supplementary Table 9f (P <= 0.05 and absolute log2FC >= 1). Main-text Dnmt3l-OE DE counts and candidate calls should be interpreted from the FDR-controlled public-count DESeq2/edgeR reanalysis. "
    "(c) Fetal human cortical plate/germinal-zone Hi-C same-TAD heatmap summarizes whether human DNMT3L-sensitive CpG-gene links fall within the same Won et al. fetal brain TAD as the candidate-gene TSS. Counts are shown for direct promoter/gene-body links, same cortical plate TAD links, same germinal-zone TAD links, and strict fetal-positive same-TAD links. This panel supports regulatory-domain plausibility for most human methylation-linked candidates but does not establish direct enhancer-promoter looping. "
    "Thus, ED Fig. 6 supports final-candidate neurodevelopmental relevance, limited overlap with the tested Dnmt3a-cKO transcriptome, and fetal cortical domain compatibility, but it does not prove methylation independence or direct methylation-expression causality. "
    "Abbreviations: OE, overexpression; cKO, conditional knockout; DE, differentially expressed; FDR, false discovery rate; logFC, log fold change; log2FC, log2 fold change; DESeq2, differential expression analysis for sequence count data 2; edgeR, empirical analysis of digital gene expression in R; GSE, Gene Expression Omnibus Series; CP, cortical plate; GZ, germinal zone; TAD, topologically associating domain; TSS, transcription start site; ED, Extended Data."
)


def replace_paragraph_text(paragraph, text: str) -> None:
    paragraph.clear()
    paragraph.add_run(text)


def main() -> None:
    doc = Document(DOCX)
    hits = {"results": False, "title": False, "legend": False}

    for para in doc.paragraphs:
        text = para.text
        if text.startswith("Because nearest-gene windows can over-assign methylation features") and "Supplementary Tables 6d, 8h, and 8i" in text:
            replace_paragraph_text(
                para,
                text.replace(
                    "(Supplementary Tables 6d, 8h, and 8i).",
                    "(Extended Data Fig. 6c; Supplementary Tables 6d, 8h, and 8i).",
                ),
            )
            hits["results"] = True
        elif text.startswith("Extended Data Figure 6. Neurogenesis annotation"):
            replace_paragraph_text(para, ED6_TITLE)
            hits["title"] = True
        elif text.startswith("Functional annotation and Dnmt3a-cKO threshold-sensitivity analyses test"):
            replace_paragraph_text(para, ED6_LEGEND)
            hits["legend"] = True

    missing = [key for key, ok in hits.items() if not ok]
    if missing:
        raise RuntimeError(f"Missing DOCX anchors: {missing}")

    doc.save(DOCX)
    print(f"Saved {DOCX}")
    print("Updated Results citation and Extended Data Figure 6 legend for panel c.")


if __name__ == "__main__":
    main()
