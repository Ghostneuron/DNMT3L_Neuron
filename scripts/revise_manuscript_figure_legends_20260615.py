#!/usr/bin/env python3
"""Revise figure legends in the 20260615 DNMT3L neuron manuscript DOCX."""

from __future__ import annotations

import shutil
from pathlib import Path

from docx import Document


MANUSCRIPT = Path("manuscript/Manuscript DNMT3L Neuron 20260615.docx")
BACKUP = Path("manuscript/Manuscript DNMT3L Neuron 20260615.before_legend_revision.docx")


PARAGRAPHS = [
    "Main Figures",
    "Main Figure 1. Target-gene expression across normal and Down syndrome brain contexts.",
    (
        "Expression and detection of DNMT3L, DNMT3A, DNMT3B, STAT1, and STAT3 were summarized across public normal "
        "developmental and Down syndrome brain datasets. (a) Human BrainSpan temporal expression shows low, episodic "
        "DNMT3L expression compared with broader DNMT3A, DNMT3B, STAT1, and STAT3 expression across prenatal and "
        "postnatal brain ages. (b) Human BrainSpan regional heatmaps summarize prenatal and adult regional variation "
        "for the same target genes. (c) Human fetal prefrontal cortex single-cell data show gestational-week trends "
        "in detection fraction and mean expression. (d) Fetal Down syndrome cortex donor-level pseudobulk effects "
        "show target-gene changes by author cell type. (e) DNMT3L detection is sparse and concentrated in selected "
        "fetal cortex cell annotations. Together, these panels support a context-specific rather than broadly "
        "upregulated DNMT3L pattern, with STAT3 showing more recurrent Down syndrome-associated signal. "
        "Abbreviations: DS, Down syndrome; RPKM, reads per kilobase of transcript per million mapped reads; TPM, "
        "transcripts per million; CPM, counts per million; pcw, post-conception weeks; mos, months; yrs, years; "
        "PFC, prefrontal cortex; DLPFC, dorsolateral PFC; VLPFC, ventrolateral PFC; ctx, cortex; SS, somatosensory; "
        "ST, superior temporal; IT, inferolateral temporal; IP, inferior parietal; MD thal, mediodorsal thalamus; "
        "AST, astrocyte; IPC, intermediate progenitor cell; IPC prol, proliferating IPC; MIC, microglia; NEU, neuron; "
        "OPC, oligodendrocyte precursor cell; RG, radial glia; VASC, vascular cell."
    ),
    "Main Figure 2. DNMT3L-linked methylation convergence nominates primary neurodevelopmental candidates.",
    (
        "Human DNMT3L overexpression, mouse Dnmt3l overexpression, Down syndrome methylation, expression, and "
        "regulatory annotation layers were integrated to prioritize candidate genes. (a) Human/mouse methylation "
        "convergence heatmap ranks genes by a heuristic integrated evidence-prioritization score combining fetal methylation-probe support, "
        "Down syndrome concordance, mouse score, and mouse-direction-concordant human expression support. (b) Top "
        "ranked evidence-prioritization candidates show the relative contribution of human methylation probes, mouse-direction-"
        "concordant human expression support, and mouse methylation score; this discovery ranking is distinct from the "
        "final 24-candidate tier order. (c) The final candidate "
        "evidence and regulatory-context heatmap summarizes methylation, expression, cross-species, and developmental "
        "regulatory support across final candidates. (d) Regulatory support counts show the number and class of "
        "overlapping human probes and mouse regulatory differentially methylated regions. (e) Regulatory enrichment "
        "tests compare primary and final candidate sets with matched background genes. The integrated evidence "
        "nominates PLPP3, FOXP1, FGFR2, ZFHX3, and SIAH3 as primary methylation-convergent candidates. Candidate "
        "nomination should be interpreted as public-data prioritization, not proof that each methylation event "
        "causally regulates the linked gene. Abbreviations: DS, Down syndrome; OE, overexpression; DMR, differentially "
        "methylated region; ATAC, assay for transposase-accessible chromatin; H3K27ac, histone H3 lysine 27 "
        "acetylation; H3K27me3, histone H3 lysine 27 trimethylation; CpG, cytosine-phosphate-guanine dinucleotide; "
        "OCR, open chromatin region."
    ),
    "Main Figure 3. Locus-level evidence for primary methylation-convergent candidates.",
    (
        "Locus panels summarize methylation, expression, and regulatory context for the five primary candidates. "
        "(a) PLPP3, (b) FOXP1, (c) FGFR2, (d) ZFHX3, and (e) SIAH3 show linked DNMT3L-sensitive methylation features "
        "and developmental regulatory annotations in human and/or mouse datasets. The mouse OE expression values shown "
        "in panel headers are gene-level public expected-count effects; for example, Zfhx3 shows only a weak "
        "non-significant mouse OE increase after transcript rows are collapsed to gene level. FOXP1 shows the clearest "
        "direction compatible with repression-like methylation-expression behavior, whereas PLPP3, FGFR2, ZFHX3, and SIAH3 "
        "show more complex regulatory methylation patterns consistent with enhancer, gene-body, timing-specific, "
        "or cell-state-dependent mechanisms. Locus proximity and chromatin annotation support regulatory plausibility "
        "but do not establish enhancer-target causality. Abbreviations: CpG, cytosine-phosphate-guanine dinucleotide; "
        "DMR, differentially methylated region; OE, overexpression; DS, Down syndrome; Dp10, Dp(10)2Yey mouse model; "
        "OCR, open chromatin region; ATAC, assay for transposase-accessible chromatin; H3K27ac, histone H3 lysine 27 "
        "acetylation; H3K4me1/3, histone H3 lysine 4 mono- or trimethylation; H3K27me3, histone H3 lysine 27 "
        "trimethylation; Mb, megabase; hg38, human genome build GRCh38; mm10, mouse genome build GRCm38."
    ),
    "Main Figure 4. Stage-specific developmental modules distinguish fetal human and postnatal mouse DNMT3L contexts.",
    (
        "Developmental module analyses were used to interpret the Dnmt3l-overexpression mouse cortex transcriptome "
        "and compare it with human Down syndrome developmental contexts. Mouse GSE268951 expression panels use public "
        "GEO processed expected counts from GSE268951_all.fpkm_mRNA_anno.txt.gz reanalyzed with DESeq2 and edgeR; "
        "FPKM is retained only as abundance context. (a) Mouse Dnmt3l-overexpression medial frontal cortex shows "
        "strongest public-count support in glial/myelination, neurodevelopmental signaling, and axon/dendrite programs. "
        "(b,c) Glial/myelination and neurodevelopmental signaling module genes are enriched for nearby Dnmt3l-sensitive "
        "methylation features. (d) Down syndrome brain concordance scores summarize whether Dnmt3l-responsive modules "
        "show similar directional behavior across mouse and human Down syndrome brain contexts. (e) Human Down syndrome "
        "neurogenesis timing modules provide fetal-stage context separate from the mouse postnatal analysis. (f) "
        "Dnmt3l-overexpression and Dnmt3a-cKO transcriptomic overlap is limited under the tested thresholds. These "
        "results support methylation-associated postnatal maturation and glial/myelination effects in mouse while "
        "preserving a separate fetal human neurodevelopment interpretation. Abbreviations: DS, Down syndrome; OE, "
        "overexpression; cKO, conditional knockout; DMR, differentially methylated region; DE, differentially expressed; "
        "FDR, false discovery rate; logFC, log fold change; FC, fold change; FPKM, fragments per kilobase of transcript "
        "per million mapped reads; P28/P39, postnatal day 28/39; mFC, medial frontal cortex; DLPFC, dorsolateral "
        "prefrontal cortex; Dp10, Dp(10)2Yey mouse model."
    ),
    "Main Figure 5. Fetal DS chromatin analyses support developmental TF landscape remodeling rather than dominant direct STAT motif enrichment.",
    (
        "Fetal Down syndrome single-cell ATAC and motif analyses were used to test whether candidate mechanisms are "
        "reflected in chromatin accessibility and transcription-factor motif landscapes. (a) Cross-context HOMER "
        "motif-family enrichment summarizes motif families across DNMT3L-overexpression, Down syndrome, and mouse "
        "contexts. (b) Ranked fetal Down syndrome peak motif shifts show distributed changes in transcription-factor "
        "motif-associated accessibility. (c) GC/content-matched chromVAR-style motif beta estimates provide a "
        "bias-controlled motif-deviation sensitivity analysis. (d) Strict peak-level differential accessibility is "
        "sparse, with significant peaks concentrated on chromosome 21 in the current analysis. (e) Core-gene cis-window "
        "accessibility summarizes local accessibility around DNMT3L, DNMT3A, DNMT3B, STAT1, STAT3, and primary "
        "candidates. Fetal Down syndrome chromatin data support altered developmental regulatory landscapes, including "
        "broader homeobox/bHLH signals and RFX/MEF2 loss, but do not support STAT as the dominant direct DNA-binding "
        "motif signal. Abbreviations: DS, Down syndrome; scATAC, single-cell assay for transposase-accessible chromatin; "
        "TF, transcription factor; HOMER, Hypergeometric Optimization of Motif EnRichment; DMR, differentially methylated "
        "region; DMC, differentially methylated cytosine or CpG; DA, differential accessibility; chromVAR, chromatin "
        "variation motif-deviation framework; PCW, post-conception week; CPM, counts per million; bHLH, basic "
        "helix-loop-helix; HMG, high-mobility group; RFX, regulatory factor X; MEF2, myocyte enhancer factor 2; STAT, "
        "signal transducer and activator of transcription; KLF/SP/GC, Kruppel-like factor/Sp-family/GC-rich motifs; "
        "CTCF/BORIS, CCCTC-binding factor/Brother of the Regulator of Imprinted Sites; ZF, zinc finger; AP1/CEBP, "
        "activator protein 1/CCAAT enhancer-binding protein; RG, radial glia; IPC, intermediate progenitor cell; EN, "
        "excitatory neuron; GABA, gamma-aminobutyric acid; OPC/OL, oligodendrocyte precursor cell/oligodendrocyte; "
        "MIC, microglia; VASC, vascular cell."
    ),
    "Main Figure 6. Integrated evidence and context limits for DNMT3L-linked mechanisms in neurodevelopment.",
    (
        "Integrated mechanism summaries place the final candidates into methylation-linked, transcription-factor/"
        "cell-state, and unresolved mechanism axes. (a) Candidate mechanism-axis heatmap summarizes evidence for "
        "DNMT3L-linked methylation, developmental transcription-factor motifs, CpG/CTCF/GC contexts, and STAT-related "
        "motif instances. (b) Candidate locus mechanism evidence summarizes locus-level support for regulatory "
        "methylation and chromatin context. (c) Fetal Down syndrome expression by mechanism axis links candidate groups "
        "to cell-type expression effects. (d) Fetal Down syndrome chromatin by mechanism axis and cell class links "
        "mechanism groups to accessibility context. Together, these data support a context-limited mechanism-evidence "
        "framework with stronger support for regulatory methylation than for direct STAT motif action. The "
        "DNMT3L-STAT1/3-HES3/ASCL1/NEUROD2/NEUROG2/CDKN1A pathway remains best supported by the direct HNP/SH-SY5Y "
        "perturbation study, not by fetal Down syndrome ATAC motif enrichment. Abbreviations: DS, Down syndrome; TF, "
        "transcription factor; CpG, cytosine-phosphate-guanine dinucleotide; CTCF, CCCTC-binding factor; GC, guanine-"
        "cytosine-rich sequence context; ATAC, assay for transposase-accessible chromatin; ctrl, control; DMR, "
        "differentially methylated region; methyl., methylation; Dev.TF, developmental transcription factor; IPC, "
        "intermediate progenitor cell; AST, astrocyte; MIC, microglia; NEU, neuron; RG, radial glia; VASC, vascular "
        "cell; EN, excitatory neuron; GABA/SST, GABAergic or somatostatin-enriched neuronal class; OPC/OL, "
        "oligodendrocyte precursor cell/oligodendrocyte."
    ),
    "Main Figure 7. Working model for DNMT3L-linked neurodevelopmental mechanisms.",
    (
        "The working model summarizes evidence boundaries for DNMT3L-linked candidate prioritization rather than a "
        "resolved causal mechanism. Public transcriptome, methylome, and chromatin datasets support candidate "
        "prioritization rather than final causal proof. Human "
        "DNMT3L-overexpression expression and methylation, mouse postnatal Dnmt3l-overexpression, Down syndrome brain "
        "expression and methylation context, and fetal chromatin/motif resources converge on primary candidates "
        "including PLPP3, FOXP1, FGFR2, ZFHX3, and SIAH3. The methylation-linked regulatory hypothesis is better "
        "supported by the current in-silico analyses, whereas STAT/TF and cell-state effects remain contextual or "
        "indirect in public Down syndrome chromatin data and require matched experimental validation. Abbreviations: "
        "DS, Down syndrome; OE, overexpression; RNA-seq, "
        "RNA sequencing; EPIC, Illumina Infinium MethylationEPIC array; WGBS, whole-genome bisulfite sequencing; cKO, "
        "conditional knockout; DMR, differentially methylated region; CpG, cytosine-phosphate-guanine dinucleotide; "
        "Dp10, Dp(10)2Yey mouse model; OCR, open chromatin region; ATAC, assay for transposase-accessible chromatin; "
        "HOMER, Hypergeometric Optimization of Motif EnRichment; chromVAR, chromatin variation motif-deviation "
        "framework; BED, Browser Extensible Data interval format; TF, transcription factor; HNP, human neural progenitor; "
        "SH-SY5Y, human neuroblastoma cell line."
    ),
    "Extended Data Figures",
    "Extended Data Figure 1. Additional mouse target-gene expression reference.",
    (
        "Mouse developmental and adult brain expression references for the five target genes are shown using Allen in "
        "situ hybridization and Allen single-cell resources. These panels complement the human BrainSpan panels in Main "
        "Figure 1 by documenting mouse developmental expression energy, adult-region enrichment, and midbrain single-cell "
        "class detection for Dnmt3l, Dnmt3a, Dnmt3b, Stat1, and Stat3. The figure supports the conservative framing that "
        "Dnmt3l is low or sparse in most normal mouse brain regions and developmental stages, whereas Dnmt3a, Dnmt3b, "
        "Stat1, and Stat3 show broader expression. Abbreviations: ISH, in situ hybridization; E, embryonic day; P, "
        "postnatal day; ABC, Allen Brain Cell atlas; MB, midbrain; GABA, gamma-aminobutyric acid; OPC, oligodendrocyte "
        "precursor cell; Oligo, oligodendrocyte; Astro-Epen, astrocyte/ependymal class; ORB, orbital area; PL, prelimbic "
        "area; AOB, accessory olfactory bulb; OB, olfactory bulb; NTS, nucleus of the solitary tract."
    ),
    "Extended Data Figure 2. Dataset-level DS expression audits and public postnatal PFC context.",
    (
        "Dataset-level Down syndrome expression audits summarize target-gene behavior across fetal, adult, organoid, "
        "induced pluripotent stem cell, neural progenitor, and public postnatal brain resources. (a,b) Target-gene "
        "effect sizes across public datasets show context-specific rather than uniform DNMT3L upregulation. (c) Fetal "
        "Down syndrome single-cell annotations provide additional cell-type context for target genes. (d) Public "
        "postnatal Down syndrome PFC Shiny visualization is included as contextual evidence because full genomic "
        "matrices require controlled access. Abbreviations: DS, Down syndrome; PFC, prefrontal cortex; DLPFC, "
        "dorsolateral PFC; HIP, hippocampus; ITC, inferior temporal cortex; OFC, orbitofrontal cortex; V1C, primary "
        "visual cortex; VFC, ventral frontal cortex; IPC, inferior parietal cortex or intermediate progenitor cell "
        "depending on panel context; MFC, medial frontal cortex; S1C, primary somatosensory cortex; STC, superior "
        "temporal cortex; iPSC, induced pluripotent stem cell; NPC, neural progenitor cell; ASO, antisense "
        "oligonucleotide; PCW, post-conception week; CPM, counts per million; TPM, transcripts per million; AST, "
        "astrocyte; RG, radial glia; OPC, oligodendrocyte precursor cell; VASC, vascular cell."
    ),
    "Extended Data Figure 3. Mouse DS model and transchromosomic expression comparators.",
    (
        "Dp10 single-cell expression, Tc1 target-kmer expression, and TcMAC21 target-coordinate audits are summarized "
        "as model-system context. These datasets test whether mouse chromosome 10 or transchromosomic Down syndrome "
        "models show target-gene changes relevant to Dnmt3l dosage. Because each model differs in genetic content, age, "
        "tissue, and assay resolution, these data are presented as comparators rather than primary evidence for the "
        "human fetal mechanism. Abbreviations: DS, Down syndrome; Dp10 or Dp(10)2Yey, mouse duplication model carrying "
        "a segment of mouse chromosome 10; WT, wild type; CPM, counts per million; log2, base-2 logarithm; Tc1 and "
        "TcMAC21, transchromosomic Down syndrome mouse models; k-mer, sequence word used for target-mapping audit."
    ),
    "Extended Data Figure 4. DNMT3L-linked methylation signatures and HOMER enrichment QC.",
    (
        "DNMT3L-overexpression methylation signatures, Down syndrome methylation concordance, and HOMER motif-enrichment "
        "outputs are summarized as quality-control and sequence-context panels. (a) Mean delta beta shows positive "
        "methylation skew in human DNMT3L-overexpression and mouse Dnmt3l-overexpression comparisons. (b) Hyper- and "
        "hypomethylated CpG/probe counts at absolute delta beta >= 0.05 show hypermethylation predominance in both "
        "overexpression datasets, especially mouse. (c) Dp10 shows higher concordance with Dnmt3l-sensitive methylation "
        "direction than Dp16 in overlapping CpGs. (d) Human convergence-filter counts summarize how DNMT3L-overexpression "
        "probes pass Down syndrome and SH-SY5Y phase-3 concordance filters. (e) Top known HOMER motifs from human "
        "DNMT3L-overexpression DMRs and adult Down syndrome brain hyper-DMRs document motif-enrichment QC and "
        "sequence-context hypotheses. Abbreviations: DS, Down syndrome; OE, overexpression; QC, quality control; CpG, "
        "cytosine-phosphate-guanine dinucleotide; DMR, differentially methylated region; DMC, differentially methylated "
        "cytosine or CpG; hyper, hypermethylated; hypo, hypomethylated; Dp10/Dp16, Dp(10)2Yey/Dp(16)Yey mouse models; "
        "HOMER, Hypergeometric Optimization of Motif EnRichment; TF, transcription factor; SH-SY5Y, human neuroblastoma "
        "cell line."
    ),
    "Extended Data Figure 5. Candidate tiering, regulatory enrichment, and mechanism-class details.",
    (
        "Expanded candidate evidence plots show final candidate regulatory classes, regulatory class counts, matched "
        "regulatory permutation results, and non-primary candidate mechanism-roadmap annotations. These panels support "
        "the main candidate prioritization while preserving the distinction between primary methylation-convergent "
        "candidates, secondary methylation/expression candidates, methylation-heavy developmental or Polycomb candidates, "
        "and expression/pathway methylation-unresolved candidates. Abbreviations: ATAC, assay for transposase-accessible "
        "chromatin; OCR, open chromatin region; DMR, differentially methylated region; H3K27ac, histone H3 lysine 27 "
        "acetylation; H3K27me3, histone H3 lysine 27 trimethylation; dev, developmental; CpG, cytosine-phosphate-"
        "guanine dinucleotide; TF, transcription factor; CTCF, CCCTC-binding factor; matched-null, matched background "
        "null set used for permutation comparison."
    ),
    "Extended Data Figure 6. Neurodevelopmental candidate annotation, Dnmt3a threshold sensitivity, and fetal Hi-C TAD support.",
    (
        "Functional annotation, threshold-sensitivity, and fetal Hi-C domain-support analyses summarize "
        "neurodevelopmental relevance among final candidates, compare public-count GSE268951 Dnmt3l-overexpression "
        "effects with the available Dnmt3a-cKO transcriptome, and test whether human DNMT3L-sensitive CpG-gene links "
        "are compatible with fetal cortical CP/GZ TADs. Sixteen of 24 final candidates are neurogenesis or "
        "neurodevelopment-related, but most are not strong single-gene differential-expression hits in every dataset. "
        "Lowering the Dnmt3a-cKO logFC threshold increases the number of Dnmt3a-responsive genes but still does not "
        "make the final candidates pass the tested transcriptomic thresholds. The Dnmt3a threshold-sensitivity count "
        "plot is a robustness visualization; its Dnmt3l-overexpression overlap line uses a nominal P <= 0.05 and "
        "absolute logFC >= 1 flag, whereas the main manuscript DE counts use FDR-controlled public-count DESeq2/edgeR "
        "results. The fetal Hi-C same-TAD panel supports regulatory-domain plausibility but not direct "
        "enhancer-promoter looping. These data support a cautious interpretation: Dnmt3l-overexpression expression "
        "effects are not simply mirrored by this Dnmt3a-cKO dataset, but methylation independence and direct "
        "methylation-expression causality are not proven. "
        "Abbreviations: OE, overexpression; cKO, conditional knockout; DE, differentially expressed; logFC, log fold "
        "change; FDR, false discovery rate; CP, cortical plate; GZ, germinal zone; TAD, topologically associating "
        "domain; TSS, transcription start site; Dnmt3l, DNA methyltransferase 3-like; Dnmt3a, DNA methyltransferase "
        "3A."
    ),
    "Extended Data Figure 7. Fetal DS ATAC and motif sensitivity analyses-I.",
    (
        "Fetal Down syndrome ATAC sensitivity panels test whether apparent motif shifts persist after accounting for "
        "peak-set composition and sequence bias. The panels include matched motif-deviation tests, GC/content-matched "
        "chromVAR-style FDR summaries, and recurrent motif peak-set comparisons for Down syndrome-open and control-open "
        "peaks. The overall pattern supports developmental transcription-factor landscape changes and RFX/MEF2 loss, "
        "while STAT motifs and specific ASCL1/NEUROG2 opening do not emerge as dominant direct signals. Abbreviations: "
        "DS, Down syndrome; ATAC, assay for transposase-accessible chromatin; FDR, false discovery rate; GC, guanine-"
        "cytosine-rich sequence context; chromVAR, chromatin variation motif-deviation framework; PCW, post-conception "
        "week; CPM, counts per million; bHLH, basic helix-loop-helix; RFX, regulatory factor X; MEF2, myocyte enhancer "
        "factor 2; STAT, signal transducer and activator of transcription."
    ),
    "Extended Data Figure 8. Fetal DS ATAC and motif sensitivity analyses-II.",
    (
        "Additional fetal Down syndrome ATAC sensitivity panels summarize motif-instance specificity, candidate local "
        "peak summaries, and mechanism-axis accessibility context. These controls distinguish broad developmental motif "
        "landscape changes from candidate-local accessibility shifts and help prevent overinterpretation of motif "
        "presence as direct transcription-factor activity. Abbreviations: DS, Down syndrome; ATAC, assay for "
        "transposase-accessible chromatin; TF, transcription factor; PCW, post-conception week; CPM, counts per million; "
        "CpG, cytosine-phosphate-guanine dinucleotide; CTCF, CCCTC-binding factor; methyl., methylation; dev, "
        "developmental; beta, regression effect estimate for the tested accessibility contrast; motif-miss, matched "
        "peaks lacking the queried motif instance."
    ),
    "Extended Data Figure 9. STAT-paper marker audit, Drosophila and early postnatal DS PFC contextual resources.",
    (
        "The STAT-paper markers HES3, ASCL1, NEUROD2, NEUROG2, and CDKN1A are audited against human DNMT3L-overexpression "
        "methylation, mouse P28 Dnmt3l-overexpression, and fetal Down syndrome chromatin analyses. (a-d) Direct gene-linked "
        "human methylation counts, mouse nearby DMR counts, mouse P28 expression effects, and binary evidence summaries "
        "separate fetal HNP methylation-independent marker evidence from later mouse methylation context. (e-h) "
        "Drosophila DNMT3L-overexpression and the Risgaard early postnatal Down syndrome PFC resource are summarized as "
        "contextual resources, including public Shiny/code availability and strict public DEG-table audit outputs. These "
        "panels document data availability and interpretive limits rather than primary mechanism evidence. Abbreviations: "
        "DS, Down syndrome; PFC, prefrontal cortex; HNP, human neural progenitor; SH-SY5Y, human neuroblastoma cell line; "
        "DMC, differentially methylated cytosine or CpG; DMR, differentially methylated region; DE, differentially expressed; "
        "DEG, differentially expressed gene; FDR, false discovery rate; OE, overexpression; cKO, conditional knockout; "
        "RNA-seq, RNA sequencing; DNMT3A/DNMT3B, DNA methyltransferase 3A/3B; dbGaP, database of Genotypes and Phenotypes."
    ),
    "Extended Data Figure 10. Drosophila DNMT3L-OE transcriptome analysis.",
    (
        "Public processed Drosophila DNMT3L-overexpression expression datasets were analyzed as methylation-independent "
        "context. (a) Count-based DESeq2 volcano plot for GSE263069 larval brain RNA-seq shows FDR-significant "
        "DNMT3L-overexpression-responsive genes from processed expected counts. (b) DESeq2 MA plot summarizes effect "
        "sizes across expression abundance. (c) Count-based module enrichment shows that broad keyword modules do not "
        "pass FDR <= 0.05 at the strict large-effect threshold, with metabolism/detoxification showing the strongest "
        "nominal enrichment. (d) Count-based DESeq2 effects for selected fly JAK/STAT and neurogenic markers show directional changes, "
        "including pros, elav, and neur, but these are not robust FDR-level marker calls. (e) The older GSE68861 "
        "whole-larva array dataset shows larger condition-level shifts in later/generational comparisons, but public "
        "normalized values are condition-level/probe-level and should remain secondary. Abbreviations: OE, overexpression; "
        "RNA-seq, RNA sequencing; DESeq2, count-based differential-expression method; MA plot, mean abundance versus "
        "log ratio plot; FDR, false discovery rate; log2FC, base-2 log fold change; FPKM, fragments per kilobase of "
        "transcript per million mapped reads; JAK/STAT, Janus kinase/signal transducer and activator of transcription; "
        "bHLH, basic helix-loop-helix."
    ),
    "Extended Data Figure 11. Cross-species DNMT3L-OE expression conservation.",
    (
        "Human fetal cortical-neuron DNMT3L-overexpression, Drosophila third-instar larval-brain DNMT3L-overexpression, "
        "and mouse P28 medial-frontal-cortex Dnmt3l-overexpression expression effects were compared at pathway-module "
        "and curated marker-family levels. The human arm uses statistics from the published GSE157080 author-supplied "
        "count-based differential-expression table, the Drosophila arm uses the count-based DESeq2 reanalysis of "
        "processed expected counts, and the mouse arm uses the public GEO processed expected-count DESeq2/edgeR "
        "reanalysis. (a) Pathway-module "
        "median log2FC values show mostly same-direction positive human-fly broad module medians after correction, with "
        "metabolism/detoxification remaining divergent. (b) Curated STAT/TF marker-family comparisons remain mixed. "
        "(c,d) Normalized nominal up/down fractions and counts summarize pathway-module directions while accounting for "
        "different module sizes and dataset gene counts. These analyses support descriptive pathway-level responsiveness "
        "and should not be interpreted as strict fly/mammal ortholog evidence or proof of conserved STAT activity. "
        "Abbreviations: OE, overexpression; P28, postnatal day 28; mFC, medial frontal cortex; DESeq2, count-based "
        "differential-expression method; log2FC, base-2 log fold change; FDR, false discovery rate; CDK, cyclin-dependent "
        "kinase; JAK/STAT, Janus kinase/signal transducer and activator of transcription; HES/Notch, hairy/enhancer of "
        "split/Notch signaling; bHLH, basic helix-loop-helix; MAPK/WNT/BMP, mitogen-activated protein kinase/Wingless-"
        "related integration site/bone morphogenetic protein pathways; TF, transcription factor."
    ),
    "",
]


def set_paragraph_text(paragraph, text: str, bold: bool = False) -> None:
    paragraph.clear()
    run = paragraph.add_run(text)
    run.bold = bold


def main() -> None:
    if not MANUSCRIPT.exists():
        raise FileNotFoundError(MANUSCRIPT)
    if not BACKUP.exists():
        shutil.copy2(MANUSCRIPT, BACKUP)
    doc = Document(MANUSCRIPT)
    if len(doc.paragraphs) < len(PARAGRAPHS):
        for _ in range(len(PARAGRAPHS) - len(doc.paragraphs)):
            doc.add_paragraph()
    for idx, text in enumerate(PARAGRAPHS):
        is_heading = text in {"Main Figures", "Extended Data Figures"} or text.startswith("Main Figure") or text.startswith(
            "Extended Data Figure"
        )
        set_paragraph_text(doc.paragraphs[idx], text, bold=is_heading and bool(text))
    for paragraph in doc.paragraphs[len(PARAGRAPHS) :]:
        set_paragraph_text(paragraph, "")
    doc.save(MANUSCRIPT)
    print(f"Revised legends in {MANUSCRIPT}")
    print(f"Backup: {BACKUP}")


if __name__ == "__main__":
    main()
