#!/usr/bin/env python3
"""Create a manuscript working-model schematic for DNMT3L candidate mechanisms."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


OUTDIR = Path("figures/manuscript_model")
NOTE = Path("notes/manuscript_model_schematic.md")


COLORS = {
    "human": "#2F6F8F",
    "mouse": "#6B8E23",
    "ds": "#B95F27",
    "chromatin": "#7B5EA7",
    "integration": "#2C3E50",
    "primary": "#D85C4A",
    "complex": "#E3A33B",
    "polycomb": "#6E5A9B",
    "validation": "#4A5568",
    "light_human": "#DCECF3",
    "light_mouse": "#E6EED9",
    "light_ds": "#F5E3D6",
    "light_chromatin": "#ECE5F6",
    "light_primary": "#F8DFDB",
    "light_complex": "#F8EAC6",
    "light_polycomb": "#E8E2F3",
    "light_validation": "#ECEFF2",
}


def add_box(
    ax,
    x: float,
    y: float,
    w: float,
    h: float,
    title: str,
    body: str,
    facecolor: str,
    edgecolor: str,
    title_color: str | None = None,
    fontsize: float = 8.6,
) -> None:
    # 1. Drop shadow (drawn first to be underneath)
    shadow = FancyBboxPatch(
        (x + 0.003, y - 0.003),
        w,
        h,
        boxstyle="round,pad=0.018,rounding_size=0.025",
        linewidth=0,
        edgecolor="none",
        facecolor="#000000",
        alpha=0.06,
        zorder=1,
    )
    ax.add_patch(shadow)

    # 2. Main card
    patch = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.018,rounding_size=0.025",
        linewidth=1.2,
        edgecolor=edgecolor,
        facecolor=facecolor,
        zorder=2,
    )
    ax.add_patch(patch)

    # 3. Separator line below title (shifted down to y + h - 0.060)
    ax.plot(
        [x + 0.018, x + w - 0.018],
        [y + h - 0.060, y + h - 0.060],
        color=edgecolor,
        linewidth=0.8,
        alpha=0.35,
        zorder=3,
    )

    # 4. Title text (aligned bottom at y + h - 0.055, using fontsize=8.2 for safe containment)
    ax.text(
        x + 0.020,
        y + h - 0.055,
        title,
        ha="left",
        va="bottom",
        fontsize=8.2,
        fontweight="bold",
        color=title_color or edgecolor,
        zorder=3,
    )

    # 5. Body text (starts at y + h - 0.075 to sit cleanly below the line)
    ax.text(
        x + 0.020,
        y + h - 0.075,
        body,
        ha="left",
        va="top",
        fontsize=fontsize,
        color="#1F2933",
        linespacing=1.18,
        zorder=3,
    )


def add_arrow(
    ax,
    start: tuple[float, float],
    end: tuple[float, float],
    color: str = "#374151",
    connectionstyle: str = "arc3,rad=0.0",
    arrowstyle: str = "-|>",
    shrinkB: float = 6,
) -> None:
    arrow = FancyArrowPatch(
        start,
        end,
        arrowstyle=arrowstyle,
        mutation_scale=11,
        linewidth=1.6,
        color=color,
        shrinkA=0,  # Starts exactly at boundary
        shrinkB=shrinkB if arrowstyle != "-" else 0,  # Shrink if there is an arrowhead
        connectionstyle=connectionstyle,
        zorder=3,
    )
    ax.add_patch(arrow)


def add_z_arrow(ax, start: tuple[float, float], end: tuple[float, float], x_mid: float, color: str = "#374151") -> None:
    from matplotlib.path import Path
    x_start, y_start = start
    x_end, y_end = end
    
    r = 0.008
    dx1 = abs(x_mid - x_start)
    dx2 = abs(x_end - x_mid)
    dy = abs(y_end - y_start)
    r = min(r, dx1 * 0.4, dx2 * 0.4, dy * 0.4)
    
    sign_x1 = 1.0 if x_mid > x_start else -1.0
    sign_y = 1.0 if y_end > y_start else -1.0
    sign_x2 = 1.0 if x_end > x_mid else -1.0
    
    verts = [
        (x_start, y_start),
        (x_mid - sign_x1 * r, y_start),
        (x_mid, y_start),
        (x_mid, y_start + sign_y * r),
        (x_mid, y_end - sign_y * r),
        (x_mid, y_end),
        (x_mid + sign_x2 * r, y_end),
        (x_end, y_end)
    ]
    codes = [
        Path.MOVETO,
        Path.LINETO,
        Path.CURVE3,
        Path.CURVE3,
        Path.LINETO,
        Path.CURVE3,
        Path.CURVE3,
        Path.LINETO
    ]
    path = Path(verts, codes)
    
    arrow = FancyArrowPatch(
        path=path,
        arrowstyle="-|>",
        mutation_scale=11,
        linewidth=1.6,
        color=color,
        shrinkA=0,
        shrinkB=6,
        zorder=3,
    )
    ax.add_patch(arrow)


def add_small_label(ax, x: float, y: float, text: str, color: str) -> None:
    ax.text(
        x,
        y,
        text,
        ha="center",
        va="center",
        fontsize=7.8,
        fontweight="bold",
        color="white",
        bbox={
            "boxstyle": "round,pad=0.22,rounding_size=0.1",
            "facecolor": color,
            "edgecolor": "none",
            "alpha": 0.95,
            "zorder": 4,
        },
        zorder=5,
    )


def make_figure() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(15.5, 8.7))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    fig.patch.set_facecolor("white")

    ax.text(
        0.5,
        0.965,
        "Working model: $\\mathit{DNMT3L}$-linked regulatory convergence across neurodevelopmental contexts",
        ha="center",
        va="top",
        fontsize=16,
        fontweight="bold",
        color="#111827",
    )
    ax.text(
        0.5,
        0.925,
        "Public transcriptome, methylome, and chromatin data support candidate prioritization, not final causal proof",
        ha="center",
        va="top",
        fontsize=10.3,
        color="#4B5563",
    )

    # Card dimensions (w_left increased to 0.265 to contain long title)
    h = 0.16
    w_left = 0.265
    w_center = 0.24
    w_right = 0.23

    # Row y-coordinates (highly symmetric 3x3 layout)
    y_top = 0.700
    y_mid = 0.460
    y_bot = 0.220

    # Column 1 (Left: Human, Mouse, Down Syndrome)
    add_box(
        ax,
        0.025,
        y_top,
        w_left,
        h,
        "Human $\\mathit{DNMT3L}$-OE expression + methylation",
        "GSE157080 RNA-seq expression\nGSE157078 EPIC methylation\nGSE168276 SH-SY5Y WGBS validation",
        COLORS["light_human"],
        COLORS["human"],
        fontsize=8.1,
    )
    add_box(
        ax,
        0.025,
        y_mid,
        w_left,
        h,
        "Mouse postnatal $\\mathit{Dnmt3l}$-OE model",
        "GSE268951 P28 RNA-seq expression\nGSE271949 P28 methylation DMRs\nGSE141587 P39 Dnmt3a-cKO reference",
        COLORS["light_mouse"],
        COLORS["mouse"],
        fontsize=8.1,
    )
    add_box(
        ax,
        0.025,
        y_bot,
        w_left,
        h,
        "Down syndrome brain/model context",
        "Expression: GSE301886, GSE59630,\nGSE305153/GSE305146; Dp10 GSE214273\nMethylation: GSE73747, GSE74486,\nGSE63347; Dp10 GSE74505",
        COLORS["light_ds"],
        COLORS["ds"],
        fontsize=8.0,
    )

    # Column 2 (Center: Methylation, Expression, Chromatin - all colored in grey)
    add_box(
        ax,
        0.37,
        y_top,
        w_center,
        h,
        "Integrated methylation features",
        "DNMT3L-sensitive CpGs/DMRs\nSH-SY5Y and DS concordance\nProximity-based CpG/DMR-gene links\nCross-species candidate ranking",
        "#EEF2F7",
        COLORS["integration"],
        fontsize=8.0,
    )
    add_box(
        ax,
        0.37,
        y_mid,
        w_center,
        h,
        "Expression and module context",
        "DNMT3L-OE and DS expression effects\nMouse postnatal maturation modules\nSTAT/TF axis: contextual, not dominant motif",
        "#EEF2F7",
        COLORS["integration"],
        fontsize=8.0,
    )
    add_box(
        ax,
        0.37,
        y_bot,
        w_center,
        h,
        "Regulatory chromatin and motifs",
        "GSE149268 fetal OCR/ATAC/histone\nGSE305146 fetal DS scATAC/HOMER/chromVAR\nGSE162170; mouse forebrain ATAC/histone refs",
        "#EEF2F7",
        COLORS["integration"],
        fontsize=7.8,
    )

    # Column 3 (Right: Candidates, Split, Other)
    add_box(
        ax,
        0.70,
        y_top,
        w_right,
        h,
        "Primary candidates",
        "$\\mathit{PLPP3}$, $\\mathit{FOXP1}$, $\\mathit{FGFR2}$,\n$\\mathit{ZFHX3}$, $\\mathit{SIAH3}$\nCross-species methylation\nand regulatory support",
        COLORS["light_primary"],
        COLORS["primary"],
        fontsize=7.8,
    )
    add_box(
        ax,
        0.70,
        y_mid,
        w_right,
        h,
        "Mechanism split",
        "$\\mathit{FOXP1}$: repression-like methylation-expression\n$\\mathit{PLPP3}$/$\\mathit{FGFR2}$/$\\mathit{ZFHX3}$/$\\mathit{SIAH3}$:\ncomplex regulatory methylation\nSTAT/TF effects: contextual/indirect",
        COLORS["light_complex"],
        COLORS["complex"],
        fontsize=7.7,
    )
    add_box(
        ax,
        0.70,
        y_bot,
        w_right,
        h,
        "Other candidate tiers",
        "Secondary methylation-expression\nDevelopmental/Polycomb targets\nExpression/pathway: unresolved,\nnot proven methylation-independent",
        COLORS["light_polycomb"],
        COLORS["polycomb"],
        fontsize=7.7,
    )

    # Bottom Validation Box (moved up to y=0.025 to prevent borderline clipping)
    add_box(
        ax,
        0.10,
        0.025,
        0.80,
        0.135,
        "Validation needed",
        "Paired methylation-expression modeling; locus-specific perturbation; enhancer-target mapping;\ntime- and cell-type-specific $\\mathit{DNMT3L}$ manipulation; rescue experiments",
        COLORS["light_validation"],
        COLORS["validation"],
        fontsize=8.0,
    )

    # ------------------- Arrow Routing -------------------

    # Visual boundaries of columns
    left_right_visual = 0.025 + w_left + 0.018    # 0.025 + 0.265 + 0.018 = 0.308
    center_left_visual = 0.37 - 0.018             # 0.37 - 0.018 = 0.352
    center_right_visual = 0.37 + w_center + 0.018 # 0.37 + 0.24 + 0.018 = 0.628
    right_left_visual = 0.70 - 0.018              # 0.70 - 0.018 = 0.682
    right_right_visual = 0.70 + w_right + 0.018   # 0.70 + 0.23 + 0.018 = 0.948

    # Row midpoints
    ym_top = y_top + h/2  # 0.780
    ym_mid = y_mid + h/2  # 0.540
    ym_bot = y_bot + h/2  # 0.300

    # 1. Left column to Center column
    # Human OE (top left) to Integrated (top center) -> Straight horizontal
    add_arrow(ax, (left_right_visual, ym_top), (center_left_visual - 0.008, ym_top), COLORS["human"])
    
    # Mouse OE (middle left) to Integrated (top center) -> Z-arrow with x_mid = 0.320
    add_z_arrow(ax, (left_right_visual, ym_mid), (center_left_visual - 0.004, ym_top - 0.05), 0.320, COLORS["mouse"])
    
    # Down syndrome (bottom left) to Expression and module context (middle center) -> Z-arrow with x_mid = 0.336
    add_z_arrow(ax, (left_right_visual, ym_bot), (center_left_visual - 0.004, ym_mid), 0.336, COLORS["ds"])

    # 2. Vertical arrows in Center column
    # Top Center -> Middle Center (arrowhead removed to be consistent with the bottom vertical arrow)
    add_arrow(ax, (0.49, y_top - 0.018), (0.49, y_mid + h + 0.018), COLORS["integration"], arrowstyle="-")
    # Middle Center -> Bottom Center (arrowhead removed)
    add_arrow(ax, (0.49, y_mid - 0.018), (0.49, y_bot + h + 0.018), COLORS["integration"], arrowstyle="-")

    # 3. Center column to Right column
    # Top Center to Top Right -> Straight horizontal
    add_arrow(ax, (center_right_visual, ym_top), (right_left_visual - 0.008, ym_top), COLORS["integration"])
    
    # Middle Center to Middle Right -> Straight horizontal
    add_arrow(ax, (center_right_visual, ym_mid), (right_left_visual - 0.008, ym_mid), COLORS["integration"])
    
    # Bottom Center to Middle Right -> Z-arrow with x_mid = 0.651, landing slightly lower
    add_z_arrow(ax, (center_right_visual, ym_bot), (right_left_visual - 0.008, ym_mid - 0.02), 0.651, COLORS["integration"])

    # 4. Vertical arrows in Right column
    # Top Right -> Middle Right
    add_arrow(ax, (0.815, y_top - 0.018), (0.815, y_mid + h + 0.018 + 0.008), COLORS["primary"])
    # Middle Right -> Bottom Right
    add_arrow(ax, (0.815, y_mid - 0.018), (0.815, y_bot + h + 0.018 + 0.008), COLORS["complex"])

    # 5. Validation arrows (routing down the right corridor)
    # Top Right -> Validation
    add_z_arrow(ax, (right_right_visual, ym_top), (0.922, 0.12), 0.958, COLORS["validation"])
    # Middle Right -> Validation
    add_z_arrow(ax, (right_right_visual, ym_mid), (0.922, 0.08), 0.968, COLORS["validation"])
    # Bottom Right -> Validation
    add_z_arrow(ax, (right_right_visual, ym_bot), (0.922, 0.04), 0.978, COLORS["validation"])
    # Bottom Center -> Validation -> Straight vertical (shrinkB=0 to prevent stray arrowhead inside bottom center box)
    add_arrow(ax, (0.49, y_bot - 0.018), (0.49, 0.025 + 0.135 + 0.018 + 0.008), COLORS["validation"], shrinkB=0)

    # 6. Label pills
    # Center vertical labels
    add_small_label(ax, 0.49, (y_top - 0.018 + y_mid + h + 0.018)/2, "expression + methylation", COLORS["integration"])
    add_small_label(ax, 0.49, (y_mid - 0.018 + y_bot + h + 0.018)/2, "BED/HOMER overlaps", COLORS["integration"])
    
    # Center-to-right labels
    add_small_label(ax, (center_right_visual + right_left_visual)/2, ym_top + 0.025, "tiering", COLORS["integration"])
    add_small_label(ax, (center_right_visual + right_left_visual)/2, ym_mid + 0.025, "context", COLORS["integration"])

    png = OUTDIR / "dnmt3l_candidate_mechanism_model.png"
    pdf = OUTDIR / "dnmt3l_candidate_mechanism_model.pdf"
    top_png = Path("figures/dnmt3l_candidate_mechanism_model.png")
    
    fig.savefig(png, dpi=260, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    fig.savefig(top_png, dpi=260, bbox_inches="tight")
    plt.close(fig)


def write_note() -> None:
    NOTE.parent.mkdir(parents=True, exist_ok=True)
    NOTE.write_text(
        """# Manuscript Model Schematic

Date: 2026-06-09

Generated files:

- `figures/manuscript_model/dnmt3l_candidate_mechanism_model.png`
- `figures/manuscript_model/dnmt3l_candidate_mechanism_model.pdf`

The schematic summarizes the current working model: DNMT3L/Dnmt3l perturbation methylation data, DS methylation/expression context, and fetal regulatory chromatin annotations nominate neurodevelopmental candidate genes. The figure intentionally frames the result as candidate prioritization and mechanism generation rather than causal proof.

Recommended use: Figure 5 or final model panel after the expression, candidate evidence, locus, and regulatory enrichment figures.
""",
        encoding="utf-8",
    )


def main() -> None:
    make_figure()
    write_note()
    print("Wrote figures/manuscript_model/dnmt3l_candidate_mechanism_model.png")
    print("Wrote figures/manuscript_model/dnmt3l_candidate_mechanism_model.pdf")
    print("Wrote notes/manuscript_model_schematic.md")


if __name__ == "__main__":
    main()
