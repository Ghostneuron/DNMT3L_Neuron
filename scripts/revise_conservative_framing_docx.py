from pathlib import Path

from docx import Document


DOCX = Path("manuscript/Manuscript DNMT3L Neuron 20260617.docx")

REPLACEMENTS = {
    "Main Figure 6. Integrated evidence for stage- and cell-type-dependent DNMT3L mechanisms in neurodevelopment.": (
        "Main Figure 6. Integrated evidence and context limits for DNMT3L-linked mechanisms in neurodevelopment."
    ),
    "Main Figure 7. Working model for DNMT3L-linked neurodevelopmental mechanisms.": (
        "Main Figure 7. Evidence map and working hypothesis for DNMT3L-linked neurodevelopmental mechanisms."
    ),
    "Extended Data Figure 9. Reconciliation of STAT-paper markers, Drosophila and early postnatal DS PFC contextual resources.": (
        "Extended Data Figure 9. STAT-paper marker audit, Drosophila, and early postnatal DS PFC contextual resources."
    ),
}

FIG6_PREFIX = (
    "Integrated mechanism summaries place the final candidates into methylation-linked, "
    "transcription-factor/cell-state, and unresolved mechanism axes."
)
FIG6_TEXT = (
    "Integrated mechanism summaries place the final candidates into methylation-linked, transcription-factor/cell-state, "
    "and unresolved mechanism axes. (a) Candidate mechanism-axis heatmap summarizes evidence for DNMT3L-linked "
    "methylation, developmental transcription-factor motifs, CpG/CTCF/GC contexts, and STAT-related motif instances. "
    "(b) Candidate locus mechanism evidence summarizes locus-level support for regulatory methylation and chromatin "
    "context. (c) Fetal Down syndrome expression by mechanism axis links candidate groups to cell-type expression "
    "effects. (d) Fetal Down syndrome chromatin by mechanism axis and cell class links mechanism groups to accessibility "
    "context. Together, these data support a context-limited mechanism-evidence framework with stronger support for "
    "regulatory methylation than for direct STAT motif action. The DNMT3L-STAT1/3-HES3/ASCL1/NEUROD2/NEUROG2/CDKN1A "
    "pathway remains best supported by the direct HNP/SH-SY5Y perturbation study, not by fetal Down syndrome ATAC motif "
    "enrichment. Abbreviations: DS, Down syndrome; TF, transcription factor; CpG, cytosine-phosphate-guanine "
    "dinucleotide; CTCF, CCCTC-binding factor; GC, guanine-cytosine-rich sequence context; ATAC, assay for "
    "transposase-accessible chromatin; ctrl, control; DMR, differentially methylated region; methyl., methylation; "
    "Dev.TF, developmental transcription factor; IPC, intermediate progenitor cell; AST, astrocyte; MIC, microglia; "
    "NEU, neuron; RG, radial glia; VASC, vascular cell; EN, excitatory neuron; GABA/SST, GABAergic or "
    "somatostatin-enriched neuronal class; OPC/OL, oligodendrocyte precursor cell/oligodendrocyte."
)

FIG7_PREFIX = "The working model summarizes the evidence boundaries for a dual-axis DNMT3L mechanism."
FIG7_TEXT = (
    "The schematic summarizes evidence boundaries for DNMT3L-linked candidate prioritization rather than a resolved "
    "causal mechanism. Public transcriptome, methylome, and chromatin datasets support candidate prioritization rather "
    "than final causal proof. Human DNMT3L-overexpression expression and methylation, mouse postnatal Dnmt3l-"
    "overexpression, Down syndrome brain expression and methylation context, and fetal chromatin/motif resources "
    "converge on primary candidates including PLPP3, FOXP1, FGFR2, ZFHX3, and SIAH3. The methylation-linked regulatory "
    "hypothesis is better supported by the current in-silico analyses, whereas STAT/TF and cell-state effects remain "
    "contextual or indirect in public Down syndrome chromatin data and require matched experimental validation. "
    "Abbreviations: DS, Down syndrome; OE, overexpression; RNA-seq, RNA sequencing; EPIC, Illumina Infinium "
    "MethylationEPIC array; WGBS, whole-genome bisulfite sequencing; cKO, conditional knockout; DMR, differentially "
    "methylated region; CpG, cytosine-phosphate-guanine dinucleotide; Dp10, Dp(10)2Yey mouse model; OCR, open "
    "chromatin region; ATAC, assay for transposase-accessible chromatin; HOMER, Hypergeometric Optimization of Motif "
    "EnRichment; chromVAR, chromatin variation motif-deviation framework; BED, Browser Extensible Data interval format; "
    "TF, transcription factor; HNP, human neural progenitor; SH-SY5Y, human neuroblastoma cell line."
)


def replace_paragraph(paragraph, text: str, *, bold=False):
    paragraph.clear()
    run = paragraph.add_run(text)
    run.bold = bold


def main():
    doc = Document(DOCX)
    for paragraph in doc.paragraphs:
        text = paragraph.text.strip()
        if text in REPLACEMENTS:
            replace_paragraph(paragraph, REPLACEMENTS[text], bold=True)
        elif text.startswith(FIG6_PREFIX):
            replace_paragraph(paragraph, FIG6_TEXT)
        elif text.startswith(FIG7_PREFIX):
            replace_paragraph(paragraph, FIG7_TEXT)
    doc.save(DOCX)


if __name__ == "__main__":
    main()
