#!/usr/bin/env python3
"""Make Qin/PubPeer-aware wording revisions in the clean submission DOCX."""

from __future__ import annotations

from pathlib import Path

from docx import Document


INFILE = Path("manuscript/Manuscript_DNMT3L_Neuron_MolNeurobio_submission_20260618_clean.docx")
OUTFILE = Path(
    "manuscript/Manuscript_DNMT3L_Neuron_MolNeurobio_submission_20260618_clean_qin_cautious.docx"
)


REPLACEMENTS = {
    (
        "However, DNMT3L biology in neurodevelopment is not limited to methylation. "
        "Qin et al. reported that DNMT3L can promote neural differentiation by enhancing STAT1 and STAT3 phosphorylation and nuclear translocation, activating downstream neurogenic factors including HES3, ASCL1, NEUROD2, NEUROG2, and CDKN1A in a manner that was not rescued by a DNA methylation inhibitor [10]. "
        "This raises an alternative or parallel model: DNMT3L may affect neural differentiation through signaling and transcription-factor activity independent of DNA methylation. "
    ): (
        "However, DNMT3L biology in neurodevelopment may not be limited to methylation. "
        "Qin et al. reported a perturbation model in which DNMT3L promoted neural differentiation with increased STAT1 and STAT3 phosphorylation/nuclear translocation and changes in downstream neurogenic factors including HES3, ASCL1, NEUROD2, NEUROG2, and CDKN1A, in a manner not rescued by a DNA methylation inhibitor [10]. "
        "This report raised an alternative or parallel model: DNMT3L may affect neural differentiation through signaling and transcription-factor activity independent of DNA methylation, although this mechanism requires independent validation in matched neurodevelopmental contexts. "
    ),
    (
        "In parallel, DNMT3L/Dnmt3l perturbation datasets in human neural progenitor-derived cells, SH-SY5Y neuronal differentiation, mouse cortex, and Drosophila provide experimental contexts in which methylation-dependent and methylation-independent signatures can be examined with explicit context limits [12,9,10]."
    ): (
        "In parallel, DNMT3L/Dnmt3l perturbation datasets in human neural progenitor-derived cells, SH-SY5Y neuronal differentiation, mouse cortex, and Drosophila provide reported experimental contexts in which methylation-dependent and methylation-independent signatures can be examined with explicit context limits [12,9,10]."
    ),
    (
        "Previous studies established important individual pieces: fetal DS methylation, DNMT3L overexpression methylation signatures, DNMT3L-STAT1/STAT3 signaling, DNMT3A-dependent neurogenic methylation, and DS single-cell multiomic dysregulation. "
    ): (
        "Previous studies established or proposed important individual pieces: fetal DS methylation, DNMT3L overexpression methylation signatures, DNMT3L-STAT1/STAT3 signaling, DNMT3A-dependent neurogenic methylation, and DS single-cell multiomic dysregulation. "
    ),
    (
        "These findings do not contradict the earlier perturbation study. Rather, they define the boundary between experimental mechanisms and public-data inference. "
        "The HNP/SH-SY5Y study directly supports a DNMT3L-STAT1/STAT3 differentiation axis in that model system. "
    ): (
        "These findings do not contradict the earlier perturbation report. Rather, they define the boundary between a model-specific perturbation mechanism and public-data inference. "
        "The HNP/SH-SY5Y study provides the primary direct support for a DNMT3L-STAT1/STAT3 differentiation axis in that system, but this axis still requires independent validation in matched fetal or DS neurodevelopmental contexts. "
    ),
    (
        "STAT/TF involvement is supported by prior direct perturbation work and by STAT3 expression relevance in DS contexts, but public DS chromatin and motif analyses do not identify STAT as the dominant direct regulatory motif signal. "
    ): (
        "STAT/TF involvement is supported mainly by a prior perturbation report and by STAT3 expression relevance in DS contexts, but public DS chromatin and motif analyses do not identify STAT as the dominant direct regulatory motif signal. "
    ),
    (
        "Prior work established fetal DS hypermethylation and implicated DNMT3L [8], demonstrated DNMT3L-induced methylation changes during SH-SY5Y neuronal differentiation [9], and showed that DNMT3L can promote neural differentiation through STAT1/STAT3 phosphorylation independent of DNA methylation in HNP/SH-SY5Y and mouse perturbation systems [10]. "
    ): (
        "Prior work established fetal DS hypermethylation and implicated DNMT3L [8], demonstrated DNMT3L-induced methylation changes during SH-SY5Y neuronal differentiation [9], and reported that DNMT3L can promote neural differentiation through STAT1/STAT3 phosphorylation independent of DNA methylation in HNP/SH-SY5Y and mouse perturbation systems [10]. "
    ),
    (
        "Qin et al. found that DNMT3L overexpression can promote neural differentiation through STAT1/STAT3 phosphorylation and downstream neurogenic markers such as HES3, ASCL1, NEUROD2, NEUROG2, and CDKN1A, with limited methylation changes in their early differentiated HNP model [10]. "
        "The present results do not contradict this mechanism. Instead, they define its boundary: the STAT1/STAT3 differentiation axis is strongly supported in the direct perturbation system, whereas public fetal DS chromatin and motif data do not show STAT motifs, ASCL1, or NEUROG2 as dominant direct chromatin-opening signals in DS fetal cortex."
    ): (
        "Qin et al. reported that DNMT3L overexpression can promote neural differentiation through STAT1/STAT3 phosphorylation and downstream neurogenic markers such as HES3, ASCL1, NEUROD2, NEUROG2, and CDKN1A, with limited methylation changes in their early differentiated HNP model [10]. "
        "The present results do not contradict this proposed mechanism. Instead, they define its current evidentiary boundary: the STAT1/STAT3 differentiation axis is best supported by the original perturbation system, whereas public fetal DS chromatin and motif data do not show STAT motifs, ASCL1, or NEUROG2 as dominant direct chromatin-opening signals in DS fetal cortex."
    ),
    (
        "The main positive evidence comes from Qin et al., who experimentally demonstrated DNMT3L interaction with STAT1/STAT3 signaling and rescue of pro-differentiation effects by STAT phosphorylation inhibitors [10]. "
    ): (
        "The main positive prior evidence comes from Qin et al., who reported DNMT3L interaction with STAT1/STAT3 signaling and rescue of pro-differentiation effects by STAT phosphorylation inhibitors [10]. "
    ),
    (
        "The Qin et al. STAT1/STAT3 differentiation axis was observed in a direct perturbation model designed to test DNMT3L function during neural differentiation, whereas the mouse Dnmt3l-OE module results emphasize postnatal glial/myelination and maturation programs. "
    ): (
        "The Qin et al. STAT1/STAT3 differentiation axis was reported in a direct perturbation model designed to test DNMT3L function during neural differentiation, whereas the mouse Dnmt3l-OE module results emphasize postnatal glial/myelination and maturation programs. "
    ),
    (
        "The evidence for methylation-independent STAT/TF mechanisms is important but currently better supported by direct perturbation studies than by public DS fetal chromatin motifs. "
    ): (
        "The evidence for methylation-independent STAT/TF mechanisms is important but currently supported mainly by prior perturbation studies rather than by public DS fetal chromatin motifs. "
    ),
    (
        "The DNMT3L-STAT1/STAT3-HES3/ASCL1/NEUROD2/NEUROG2/CDKN1A pathway remains best supported by the direct HNP/SH-SY5Y perturbation study, not by fetal Down syndrome ATAC motif enrichment. "
    ): (
        "The DNMT3L-STAT1/STAT3-HES3/ASCL1/NEUROD2/NEUROG2/CDKN1A pathway remains supported chiefly by the HNP/SH-SY5Y perturbation report and was not independently recovered as a dominant fetal Down syndrome ATAC motif-enrichment signal. "
    ),
}


def replace_preserving_first_run(paragraph, old: str, new: str) -> bool:
    text = paragraph.text
    if old not in text:
        return False
    # Most manuscript paragraphs are a single formatted run after citation-field unlinking.
    # Preserve paragraph-level style and the first run's character formatting.
    runs = paragraph.runs
    fmt = None
    if runs:
        first = runs[0]
        fmt = {
            "bold": first.bold,
            "italic": first.italic,
            "underline": first.underline,
            "font_name": first.font.name,
            "font_size": first.font.size,
        }
    new_text = text.replace(old, new)
    for run in list(runs):
        run.text = ""
    run = paragraph.add_run(new_text)
    if fmt:
        run.bold = fmt["bold"]
        run.italic = fmt["italic"]
        run.underline = fmt["underline"]
        run.font.name = fmt["font_name"]
        run.font.size = fmt["font_size"]
    return True


def main() -> None:
    doc = Document(INFILE)
    counts = {old: 0 for old in REPLACEMENTS}
    for paragraph in doc.paragraphs:
        for old, new in REPLACEMENTS.items():
            if replace_preserving_first_run(paragraph, old, new):
                counts[old] += 1
    missed = [old for old, count in counts.items() if count != 1]
    if missed:
        raise SystemExit(f"Unexpected replacement counts: {[(old[:80], counts[old]) for old in missed]}")
    doc.save(OUTFILE)
    print(f"Wrote {OUTFILE}")


if __name__ == "__main__":
    main()
