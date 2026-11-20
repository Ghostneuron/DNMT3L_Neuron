"""Render the two revised figures as vector PDFs from results.json."""
import json
from pathlib import Path
import argparse

from reportlab.pdfgen import canvas
from reportlab.lib.colors import HexColor, white

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("results", type=Path, help="results.json written by reanalyze_r1.py")
parser.add_argument("--outdir", type=Path, default=Path("r1_figures"))
args = parser.parse_args()
args.outdir.mkdir(parents=True, exist_ok=True)
data = json.loads(args.results.read_text())
ink = HexColor("#17283A")
muted = HexColor("#485969")
grid = HexColor("#D7DEE4")
blue = HexColor("#18697C")
grey = HexColor("#778693")


def label(c, x, y, message, size=10, bold=False, color=ink, align="left"):
    c.setFont("Helvetica-Bold" if bold else "Helvetica", size)
    c.setFillColor(color)
    command = {"left": c.drawString, "center": c.drawCentredString, "right": c.drawRightString}
    command[align](x, y, message)


def line(c, x, y, xx, yy, color=grid, width=1):
    c.setStrokeColor(color)
    c.setLineWidth(width)
    c.line(x, y, xx, yy)


def circle(c, x, y, col, filled):
    c.setStrokeColor(col)
    c.setFillColor(col if filled else white)
    c.setLineWidth(1.2)
    c.circle(x, y, 3.6, stroke=1, fill=1)


def save(c):
    c.showPage()
    c.save()


w, h = 525, 380
c = canvas.Canvas(str(args.outdir / "Figure_1_Concordance.pdf"), pagesize=(w, h), invariant=1)
label(c, 14, 356, "A   Human methylation concordance", 12, True)
circle(c, 18, 333, blue, True)
label(c, 28, 330, "Observed")
circle(c, 98, 333, grey, False)
label(c, 108, 330, "Expected within comparable CpGs")
label(c, 506, 330, "Excess (pp)", 9, True, align="right")
contexts = ("GSE73747  Fetal, pooled", "GSE73747  Fetal frontal",
            "GSE73747  Fetal temporal", "GSE74486  Fetal cerebrum",
            "GSE74486  Adult sorted neurons", "GSE74486  Adult sorted glia")
left, right = 206, 428
def x_fraction(value):
    return left + (value - .45) / .25 * (right - left)
for value in (.45, .50, .55, .60, .65, .70):
    px = x_fraction(value)
    line(c, px, 310, px, 180)
    label(c, px, 164, f"{value:.0%}", 9, color=muted, align="center")
for index, (name, row) in enumerate(zip(contexts, data["human"])):
    y = 303 - 22 * index
    obs, exp = row["observed_agreement"], row["expected_stratified"]
    label(c, 14, y - 3, name, 9.5)
    line(c, x_fraction(obs), y, x_fraction(exp), y, grey, 1.6)
    circle(c, x_fraction(exp), y, grey, False)
    circle(c, x_fraction(obs), y, blue, True)
    label(c, 506, y - 3, f'{row["excess_stratified_pp"]:+.2f}', 9.5, align="right")
line(c, 14, 143, 506, 143)
label(c, 14, 122, "B   Mouse segmental-trisomy comparison", 12, True)
label(c, 14, 99, "CpG set", 9, True, muted)
label(c, 340, 99, "Dp(10)1Yey", 9, True, muted, "right")
label(c, 423, 99, "Dp(16)1Yey", 9, True, muted, "right")
label(c, 506, 99, "Gap (pp)", 9, True, muted, "right")
pair = data["mouse"][0]
original = {r["model"]: 100 * float(r["concordant_fraction"]) for r in data["mouse_unpaired"]}
for y, title, a, b in [(76, "Model-specific overlaps", original["Dp10"], original["Dp16"]),
                        (52, "Same 3,089 CpGs in both models", 100 * pair["dp10_agreement"], 100 * pair["dp16_agreement"])]:
    label(c, 14, y, title)
    label(c, 340, y, f"{a:.1f}%", 10, align="right")
    label(c, 423, y, f"{b:.1f}%", 10, align="right")
    label(c, 506, y, f"{a-b:+.1f}", 10, True, blue, "right")
line(c, 14, 66, 506, 66)
label(c, 14, 23, "Agreement is directional. CpGs are genomic loci. pp, percentage points.", 9, color=muted)
save(c)

c = canvas.Canvas(str(args.outdir / "Figure_2_Rank_sensitivity.pdf"), pagesize=(525, 314), invariant=1)
label(c, 14, 290, "Overlap with the corrected ranking", 12, True)
labels = {
    "submitted_st4c_score": "Submitted discovery score",
    "equal_weight": "Equal weight per component",
    "no_expression": "Methylation components only",
    "no_human_methylation": "Remove human methylation",
    "no_mouse_methylation": "Remove mouse methylation",
    "no_dp10": "Remove Dp(10)2Yey expression",
    "no_distal_mouse": "Remove distal mouse DMR points",
}
left, right = 220, 479
for index in range(0, 25, 6):
    x = left + index / 24 * (right - left)
    line(c, x, 270, x, 62)
    label(c, x, 48, str(index), 9, color=muted, align="center")
rows = [row for row in data["ranks"]["sensitivity"] if row["scenario"] in labels]
for i, row in enumerate(rows):
    y = 256 - 29 * i
    n = row["top24_overlap_with_corrected"]
    label(c, 14, y-3, labels[row["scenario"]], 9.8)
    c.setFillColor(blue if i > 0 else grey)
    x = left + n / 24 * (right - left)
    c.rect(left, y-6, x-left, 12, stroke=0, fill=1)
    label(c, x+7, y-3, str(n), 9, True)
label(c, 348, 24, "Shared genes among the 24 highest-scoring gene IDs", 9, color=muted, align="center")
save(c)
