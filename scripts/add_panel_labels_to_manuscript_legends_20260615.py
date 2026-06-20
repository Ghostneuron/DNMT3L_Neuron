#!/usr/bin/env python3
"""Add explicit panel-by-panel labels to revised DNMT3L manuscript legends."""

from __future__ import annotations

import shutil
from pathlib import Path

from docx import Document


MANUSCRIPT = Path("manuscript/Manuscript DNMT3L Neuron 20260615.docx")
BACKUP = Path("manuscript/Manuscript DNMT3L Neuron 20260615.before_panel_label_revision.docx")


REPLACEMENTS = {
    "Main Figure 4. Stage-specific developmental modules distinguish fetal human and postnatal mouse DNMT3L contexts.": (
        "Developmental module analyses were used to interpret the Dnmt3l-overexpression mouse cortex transcriptome "
        "and compare it with human Down syndrome developmental contexts. Mouse GSE268951 expression panels use the "
        "public GEO processed expected counts from GSE268951_all.fpkm_mRNA_anno.txt.gz reanalyzed with DESeq2 and "
        "edgeR; FPKM is retained only as abundance context. (a) Mouse Dnmt3l-overexpression medial frontal cortex "
        "shows strongest public-count support in glial/myelination, neurodevelopmental signaling, and axon/dendrite programs. "
        "(b) Glial/myelination module genes are enriched for nearby Dnmt3l-sensitive methylation features. "
        "(c) Neurodevelopmental signaling module genes are also enriched for nearby Dnmt3l-sensitive methylation "
        "features. (d) Down syndrome brain concordance scores summarize whether Dnmt3l-responsive modules show similar "
        "directional behavior across mouse and human Down syndrome brain contexts. (e) Human Down syndrome "
        "neurogenesis timing modules provide fetal-stage context separate from the mouse postnatal analysis. "
        "(f) Dnmt3l-overexpression and Dnmt3a-cKO transcriptomic overlap is limited under the tested thresholds. "
        "These results support methylation-associated postnatal maturation and glial/myelination effects in mouse while "
        "preserving a separate fetal human neurodevelopment interpretation. Abbreviations: DS, Down syndrome; OE, "
        "overexpression; cKO, conditional knockout; DMR, differentially methylated region; DE, differentially expressed; "
        "FDR, false discovery rate; logFC, log fold change; FC, fold change; FPKM, fragments per kilobase of transcript "
        "per million mapped reads; P28/P39, postnatal day 28/39; mFC, medial frontal cortex; DLPFC, dorsolateral "
        "prefrontal cortex; Dp10, Dp(10)2Yey mouse model."
    ),
    "Extended Data Figure 1. Additional mouse target-gene expression reference.": (
        "Mouse developmental and adult brain expression references for the five target genes are shown using Allen "
        "Institute resources. (a) Developing mouse brain in situ hybridization expression-energy trajectories summarize "
        "Dnmt3l, Dnmt3a, Dnmt3b, Stat1, and Stat3 across embryonic and postnatal ages. (b) Adult mouse in situ "
        "hybridization top-region plots rank the regions with strongest expression energy for each target gene. "
        "(c) Allen single-cell midbrain detection fractions show target-gene detection across broad midbrain cell "
        "classes. These panels complement the human BrainSpan panels in Main Figure 1 and support the conservative "
        "framing that Dnmt3l is low or sparse in most normal mouse brain regions and developmental stages, whereas "
        "Dnmt3a, Dnmt3b, Stat1, and Stat3 show broader expression. Abbreviations: ISH, in situ hybridization; E, "
        "embryonic day; P, postnatal day; ABC, Allen Brain Cell atlas; MB, midbrain; GABA, gamma-aminobutyric acid; "
        "OPC, oligodendrocyte precursor cell; Oligo, oligodendrocyte; Astro-Epen, astrocyte/ependymal class; ORB, "
        "orbital area; PL, prelimbic area; AOB, accessory olfactory bulb; OB, olfactory bulb; NTS, nucleus of the "
        "solitary tract."
    ),
    "Extended Data Figure 2. Dataset-level DS expression audits and public postnatal PFC context.": (
        "Dataset-level Down syndrome expression audits summarize target-gene behavior across fetal, adult, organoid, "
        "induced pluripotent stem cell, neural progenitor, and public postnatal brain resources. (a) Target-gene effect "
        "sizes across public datasets show context-specific rather than uniform DNMT3L upregulation. (b) Ranked "
        "dataset-context effects summarize the strongest positive and negative DNMT3L-associated signals across the "
        "surveyed resources. (c) Fetal Down syndrome single-cell subclass effects provide additional cell-type context "
        "for DNMT3L and related target-gene signals. (d) Public postnatal Down syndrome PFC Shiny visualization is "
        "included as contextual evidence because full genomic matrices require controlled access. Abbreviations: DS, "
        "Down syndrome; PFC, prefrontal cortex; DLPFC, dorsolateral PFC; HIP, hippocampus; ITC, inferior temporal "
        "cortex; OFC, orbitofrontal cortex; V1C, primary visual cortex; VFC, ventral frontal cortex; IPC, inferior "
        "parietal cortex or intermediate progenitor cell depending on panel context; MFC, medial frontal cortex; S1C, "
        "primary somatosensory cortex; STC, superior temporal cortex; iPSC, induced pluripotent stem cell; NPC, neural "
        "progenitor cell; ASO, antisense oligonucleotide; PCW, post-conception week; CPM, counts per million; TPM, "
        "transcripts per million; AST, astrocyte; RG, radial glia; OPC, oligodendrocyte precursor cell; VASC, vascular "
        "cell."
    ),
    "Extended Data Figure 3. Mouse DS model and transchromosomic expression comparators.": (
        "Mouse chromosome-10 duplication and transchromosomic Down syndrome model resources are summarized as "
        "model-system context. (a) Dp10 single-cell pseudobulk target-gene log2(CPM+1) values compare wild-type and "
        "Dp(10)2Yey samples. (b) Dp(10)2Yey minus wild-type delta log2(CPM+1) summarizes target-gene direction. "
        "(c) Dp10 target-gene detection percentages show whether apparent expression changes reflect detectable cells. "
        "(d) Tc1 and TcMAC21 target-kmer CPM audits summarize human and mouse target-coordinate signal in "
        "transchromosomic resources. Because each model differs in genetic content, age, tissue, and assay resolution, "
        "these data are presented as comparators rather than primary evidence for the human fetal mechanism. "
        "Abbreviations: DS, Down syndrome; Dp10 or Dp(10)2Yey, mouse duplication model carrying a segment of mouse "
        "chromosome 10; WT, wild type; CPM, counts per million; log2, base-2 logarithm; Tc1 and TcMAC21, "
        "transchromosomic Down syndrome mouse models; k-mer, sequence word used for target-mapping audit."
    ),
    "Extended Data Figure 5. Candidate tiering, regulatory enrichment, and mechanism-class details.": (
        "Expanded candidate evidence plots support the main candidate prioritization. (a) Regulatory-context class "
        "counts summarize how many final candidates have human, mouse, or combined regulatory support. (b) Mechanism "
        "class counts summarize the distribution of methylation-expression, developmental/Polycomb, secondary, and "
        "unresolved candidate classes. (c) Matched regulatory permutation tests compare observed candidate regulatory "
        "overlap with matched-null background expectations. (d) Non-primary candidate mechanism-roadmap annotations "
        "separate repression/derepression, enhancer/gene-body, Polycomb/H3K27me3, developmental timing, pathway/network, "
        "cell-type-specificity, cross-species, and regulatory-chromatin support. These panels preserve the distinction "
        "between primary methylation-convergent candidates, secondary methylation/expression candidates, methylation-heavy "
        "developmental or Polycomb candidates, and expression/pathway methylation-unresolved candidates. Abbreviations: "
        "ATAC, assay for transposase-accessible chromatin; OCR, open chromatin region; DMR, differentially methylated "
        "region; H3K27ac, histone H3 lysine 27 acetylation; H3K27me3, histone H3 lysine 27 trimethylation; dev, "
        "developmental; CpG, cytosine-phosphate-guanine dinucleotide; TF, transcription factor; CTCF, CCCTC-binding "
        "factor; matched-null, matched background null set used for permutation comparison."
    ),
    "Extended Data Figure 6. Neurodevelopmental candidate annotation, Dnmt3a threshold sensitivity, and fetal Hi-C TAD support.": (
        "Functional annotation, threshold-sensitivity, and fetal Hi-C domain-support analyses summarize "
        "neurodevelopmental relevance among final candidates, compare public-count GSE268951 Dnmt3l-overexpression "
        "effects with the available Dnmt3a-cKO transcriptome, and test whether human DNMT3L-sensitive CpG-gene links "
        "are compatible with fetal cortical 3D domains. (a) Final candidate neurogenesis annotation heatmap summarizes "
        "core neurogenesis, brain-region patterning, neurite/synapse/migration, glial context, neurogenesis signaling, "
        "Dnmt3l-overexpression differential expression, and Dnmt3a-cKO differential expression annotations. (b) "
        "Dnmt3a-cKO threshold-sensitivity curves show how protein-coding Dnmt3a-responsive genes, neurogenesis-related "
        "genes, and overlap with Dnmt3l-OE differentially expressed genes change as the absolute logFC threshold is "
        "relaxed. (c) Fetal CP/GZ Hi-C same-TAD heatmap summarizes direct promoter/gene-body links, same cortical "
        "plate TAD links, same germinal-zone TAD links, and strict fetal-positive same-TAD links for human "
        "DNMT3L-sensitive candidate CpG-gene records. Sixteen of 24 final candidates are neurogenesis or "
        "neurodevelopment-related, but most are not strong single-gene differential-expression hits in every dataset. "
        "These data support a cautious interpretation: Dnmt3l-overexpression expression effects are not simply mirrored "
        "by this Dnmt3a-cKO dataset, and fetal Hi-C TAD support strengthens regulatory-domain plausibility but does not "
        "prove methylation independence, direct enhancer-promoter looping, or methylation-expression causality. "
        "Abbreviations: OE, overexpression; cKO, conditional knockout; DE, differentially expressed; logFC, log fold "
        "change; FDR, false discovery rate; CP, cortical plate; GZ, germinal zone; TAD, topologically associating "
        "domain; TSS, transcription start site; Dnmt3l, DNA methyltransferase 3-like; Dnmt3a, DNA methyltransferase "
        "3A."
    ),
    "Extended Data Figure 7. Fetal DS ATAC and motif sensitivity analyses-I.": (
        "Fetal Down syndrome ATAC sensitivity panels test whether apparent motif shifts persist after accounting for "
        "peak-set composition and sequence bias. (a) Matched motif-deviation heatmap summarizes DS-control delta "
        "motif/matched log2 ratios across recurrent motif families and fetal cortex cell classes. (b) GC/content-matched "
        "chromVAR-style signed FDR heatmap tests whether motif-deviation signals remain significant after matching. "
        "(c) Recurrent motif peak-set accessibility heatmap summarizes DS-control delta log2(peak CPM+1) for DS-open "
        "and control-open recurrent motif sets. The overall pattern supports developmental transcription-factor "
        "landscape changes and RFX/MEF2 loss, while STAT motifs and specific ASCL1/NEUROG2 opening do not emerge as "
        "dominant direct signals. Abbreviations: DS, Down syndrome; ATAC, assay for transposase-accessible chromatin; "
        "FDR, false discovery rate; GC, guanine-cytosine-rich sequence context; chromVAR, chromatin variation motif-"
        "deviation framework; PCW, post-conception week; CPM, counts per million; bHLH, basic helix-loop-helix; RFX, "
        "regulatory factor X; MEF2, myocyte enhancer factor 2; STAT, signal transducer and activator of transcription."
    ),
    "Extended Data Figure 8. Fetal DS ATAC and motif sensitivity analyses-II.": (
        "Additional fetal Down syndrome ATAC sensitivity panels distinguish broad developmental motif landscape changes "
        "from candidate-local accessibility shifts. (a) Motif-instance specificity heatmap summarizes hit delta minus "
        "motif-miss delta for recurrent motif sets across fetal cortex cell classes. (b) Candidate local peak heatmap "
        "summarizes mean PCW-adjusted Down syndrome beta across local peaks near target and candidate genes. "
        "(c) Mechanism-axis local peak heatmap summarizes mean local peak beta for candidate mechanism classes across "
        "cell classes. These controls help prevent overinterpretation of motif presence as direct transcription-factor "
        "activity. Abbreviations: DS, Down syndrome; ATAC, assay for transposase-accessible chromatin; TF, transcription "
        "factor; PCW, post-conception week; CPM, counts per million; CpG, cytosine-phosphate-guanine dinucleotide; "
        "CTCF, CCCTC-binding factor; methyl., methylation; dev, developmental; beta, regression effect estimate for the "
        "tested accessibility contrast; motif-miss, matched peaks lacking the queried motif instance."
    ),
    "Extended Data Figure 9. Reconciliation of STAT-paper markers, Drosophila and early postnatal DS PFC contextual resources.": (
        "The STAT-paper markers HES3, ASCL1, NEUROD2, NEUROG2, and CDKN1A are audited against human DNMT3L-overexpression "
        "methylation, mouse P28 Dnmt3l-overexpression, and fetal Down syndrome chromatin analyses, with non-mammalian "
        "and postnatal context shown separately. (a) Direct gene-linked human DMC/DMR counts and mouse promoter/gene-body "
        "DMR counts summarize marker methylation support. (b) Nearby mouse DMR counts within 100 kb summarize later-stage "
        "methylation context around the same markers. (c) Mouse P28 Dnmt3l-overexpression log2FC values show that marker "
        "expression effects do not pass FDR significance. (d) Binary evidence summary separates absent direct human "
        "candidate DMC/DMR support from present mouse nearby DMR support and absent FDR-significant mouse Dnmt3l-OE or "
        "Dnmt3a-cKO expression effects. (e) Drosophila dataset scoring compares GSE263069 larval brain RNA-seq and "
        "GSE68861 whole-larva array resources. (f) Risgaard early postnatal Down syndrome PFC resource availability "
        "summarizes public metadata, controlled matrices, code, Shiny visualization, and strict Shiny DEG audit outputs. "
        "(g) Qualitative public Shiny target-gene detection scores summarize DNMT3L, DNMT3A, DNMT3B, STAT1, and STAT3. "
        "(h) Strict public Shiny DEG-table queries show no target/candidate matches under default filters. These panels "
        "document data availability and interpretive limits rather than primary mechanism evidence. Abbreviations: DS, "
        "Down syndrome; PFC, prefrontal cortex; HNP, human neural progenitor; SH-SY5Y, human neuroblastoma cell line; "
        "DMC, differentially methylated cytosine or CpG; DMR, differentially methylated region; DE, differentially "
        "expressed; DEG, differentially expressed gene; FDR, false discovery rate; OE, overexpression; cKO, conditional "
        "knockout; RNA-seq, RNA sequencing; DNMT3A/DNMT3B, DNA methyltransferase 3A/3B; dbGaP, database of Genotypes "
        "and Phenotypes."
    ),
    "Extended Data Figure 11. Cross-species DNMT3L-OE expression conservation.": (
        "Human fetal cortical-neuron DNMT3L-overexpression, Drosophila third-instar larval-brain DNMT3L-overexpression, "
        "and mouse P28 medial-frontal-cortex Dnmt3l-overexpression expression effects were compared at pathway-module "
        "and curated marker-family levels. The human arm uses statistics from the published GSE157080 author-supplied "
        "count-based differential-expression table, the Drosophila arm uses the count-based DESeq2 reanalysis of "
        "processed expected counts, and the "
        "mouse arm uses the public GEO processed expected-count DESeq2/edgeR reanalysis. (a) Pathway-module "
        "median log2FC values show mostly same-direction positive human-fly broad module medians after correction, with "
        "metabolism/detoxification remaining divergent. (b) Curated STAT/TF marker-family median log2FC comparisons "
        "remain mixed across species and developmental stages. (c) Normalized nominal up/down fractions summarize "
        "pathway-module directions while accounting for different module sizes. (d) Nominal up/down gene counts show "
        "the underlying gene-number contribution to the directional fraction summaries. These analyses support "
        "descriptive pathway-level responsiveness and should not be interpreted as strict fly/mammal ortholog evidence "
        "or proof of conserved STAT activity. Abbreviations: OE, overexpression; P28, postnatal day 28; mFC, medial "
        "frontal cortex; DESeq2, count-based differential-expression method; log2FC, base-2 log fold change; FDR, false "
        "discovery rate; CDK, cyclin-dependent kinase; JAK/STAT, Janus kinase/signal transducer and activator of "
        "transcription; HES/Notch, hairy/enhancer of split/Notch signaling; bHLH, basic helix-loop-helix; MAPK/WNT/BMP, "
        "mitogen-activated protein kinase/Wingless-related integration site/bone morphogenetic protein pathways; TF, "
        "transcription factor."
    ),
}


def replace_text(paragraph, text: str) -> None:
    paragraph.clear()
    paragraph.add_run(text)


def main() -> None:
    if not BACKUP.exists():
        shutil.copy2(MANUSCRIPT, BACKUP)
    doc = Document(MANUSCRIPT)
    replacements_done = 0
    for idx, paragraph in enumerate(doc.paragraphs[:-1]):
        heading = paragraph.text.strip()
        if heading in REPLACEMENTS:
            replace_text(doc.paragraphs[idx + 1], REPLACEMENTS[heading])
            replacements_done += 1
    doc.save(MANUSCRIPT)
    print(f"Updated {replacements_done} legends in {MANUSCRIPT}")
    print(f"Backup: {BACKUP}")


if __name__ == "__main__":
    main()
