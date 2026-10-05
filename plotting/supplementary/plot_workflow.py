#!/usr/bin/env python3
"""Create an editable workflow diagram for the SoilSynCom GEM manuscript."""

from __future__ import annotations

import importlib.util
import math
import sys
import textwrap
from pathlib import Path

import matplotlib

if importlib.util.find_spec("cairo") or importlib.util.find_spec("cairocffi"):
    matplotlib.use("cairo")
else:
    matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib import patches
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from publication_style import apply_publication_style  # noqa: E402


OUT_DIR = ROOT / "outputs" / "workflow"
OUT_STEM = OUT_DIR / "S04_workflow"

PHYLUM_COLORS = {
    "Actinobacteria": "#1b9e77",
    "Alphaproteobacteria": "#d95f02",
    "Betaproteobacteria": "#7570b3",
    "Gammaproteobacteria": "#e7298a",
    "Firmicutes": "#66a61e",
    "Bacteroidetes": "#e6ab02",
}

TEXT = "#1F1F1F"
MUTED = "#6B6B6B"
EDGE = "#303030"
LIGHT_EDGE = "#9A9A9A"
INPUT_FILL = "#F4F7F8"
MODEL_FILL = "#F6F2EC"
ANALYSIS_FILL = "#F2F6F3"
OUTPUT_FILL = "#F7F4F8"
ENV_FILL = "#F2F6FA"
ARROW = "#575757"


def mm_to_in(width_mm: float, height_mm: float) -> tuple[float, float]:
    return width_mm / 25.4, height_mm / 25.4


def add_text(
    ax,
    x: float,
    y: float,
    text: str,
    *,
    size: float = 8,
    weight: str = "normal",
    ha: str = "center",
    va: str = "center",
    color: str = TEXT,
    linespacing: float = 1.15,
) -> None:
    ax.text(
        x,
        y,
        text,
        ha=ha,
        va=va,
        fontsize=size,
        fontweight=weight,
        color=color,
        linespacing=linespacing,
    )


def rounded_box(
    ax,
    x: float,
    y: float,
    w: float,
    h: float,
    *,
    facecolor: str,
    edgecolor: str = EDGE,
    linewidth: float = 0.75,
    radius: float = 0.6,
    zorder: int = 1,
) -> FancyBboxPatch:
    patch = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle=f"round,pad=0.24,rounding_size={radius}",
        facecolor=facecolor,
        edgecolor=edgecolor,
        linewidth=linewidth,
        zorder=zorder,
    )
    ax.add_patch(patch)
    return patch


def module_box(
    ax,
    x: float,
    y: float,
    w: float,
    h: float,
    title: str,
    body: str,
    *,
    facecolor: str,
    title_size: float = 8.1,
    body_size: float = 7.1,
    title_y: float | None = None,
    body_y: float | None = None,
) -> None:
    rounded_box(ax, x, y, w, h, facecolor=facecolor, edgecolor=LIGHT_EDGE, linewidth=0.65)
    if title_y is None:
        title_y = y + h - 2.1
    if body_y is None:
        body_y = y + h * 0.45
    add_text(ax, x + w / 2, title_y, title, size=title_size, weight="bold")
    add_text(ax, x + w / 2, body_y, body, size=body_size, color=TEXT)


def arrow(
    ax,
    start: tuple[float, float],
    end: tuple[float, float],
    *,
    rad: float = 0.0,
    color: str = ARROW,
    lw: float = 0.9,
    mutation_scale: float = 8,
    zorder: int = 3,
) -> None:
    patch = FancyArrowPatch(
        start,
        end,
        arrowstyle="-|>",
        connectionstyle=f"arc3,rad={rad}",
        linewidth=lw,
        color=color,
        mutation_scale=mutation_scale,
        shrinkA=1.5,
        shrinkB=1.5,
        zorder=zorder,
    )
    ax.add_patch(patch)


def mini_arrow(ax, x1: float, y1: float, x2: float, y2: float) -> None:
    arrow(ax, (x1, y1), (x2, y2), lw=0.7, mutation_scale=6)


def draw_syncom_icon(ax, x: float, y: float, w: float, h: float) -> None:
    colors = (
        [PHYLUM_COLORS["Alphaproteobacteria"]] * 6
        + [PHYLUM_COLORS["Gammaproteobacteria"]] * 5
        + [PHYLUM_COLORS["Actinobacteria"]] * 5
        + [PHYLUM_COLORS["Bacteroidetes"]] * 3
        + [PHYLUM_COLORS["Betaproteobacteria"]] * 2
        + [PHYLUM_COLORS["Firmicutes"]] * 1
    )
    rows, cols = 4, 6
    dx = w / (cols + 0.35)
    dy = h / (rows + 0.25)
    r = min(dx, dy) * 0.32
    for i, color in enumerate(colors[:21]):
        row = i // cols
        col = i % cols
        cx = x + dx * (col + 0.72)
        cy = y + h - dy * (row + 0.65)
        circle = patches.Circle((cx, cy), r, facecolor=color, edgecolor="white", linewidth=0.55)
        ax.add_patch(circle)
    rounded_box(ax, x + w * 0.08, y + h * 0.02, w * 0.84, h * 0.14, facecolor="#FFFFFF", edgecolor="#D0D0D0", linewidth=0.4, radius=0.25)
    add_text(ax, x + w / 2, y + h * 0.09, "21 strains", size=7.1, weight="bold")


def draw_genome_icon(ax, x: float, y: float, w: float, h: float) -> None:
    xs = [x + w * (0.15 + 0.7 * i / 40) for i in range(41)]
    mid = y + h * 0.52
    amp = h * 0.22
    wave1 = [(xx, mid + amp * math.sin((i / 40) * 2.0 * math.pi)) for i, xx in enumerate(xs)]
    wave2 = [(xx, mid - amp * math.sin((i / 40) * 2.0 * math.pi)) for i, xx in enumerate(xs)]
    ax.add_line(Line2D([p[0] for p in wave1], [p[1] for p in wave1], color="#4B7FAD", linewidth=1.0))
    ax.add_line(Line2D([p[0] for p in wave2], [p[1] for p in wave2], color="#4B7FAD", linewidth=1.0))
    for i in range(2, 40, 5):
        ax.add_line(Line2D([wave1[i][0], wave2[i][0]], [wave1[i][1], wave2[i][1]], color="#8AA9C5", linewidth=0.7))


def draw_plate_icon(ax, x: float, y: float, w: float, h: float) -> None:
    rounded_box(ax, x + w * 0.16, y + h * 0.18, w * 0.68, h * 0.62, facecolor="#FFFFFF", edgecolor="#8A8A8A", linewidth=0.55, radius=0.25)
    for row in range(3):
        for col in range(4):
            cx = x + w * (0.27 + col * 0.15)
            cy = y + h * (0.30 + row * 0.16)
            color = "#9E4A4F" if (row + col) % 3 == 0 else "#D7E5E6"
            ax.add_patch(patches.Circle((cx, cy), w * 0.035, facecolor=color, edgecolor="#666666", linewidth=0.35))


def draw_metabolomics_icon(ax, x: float, y: float, w: float, h: float) -> None:
    base_y = y + h * 0.25
    ax.add_line(Line2D([x + w * 0.16, x + w * 0.84], [base_y, base_y], color="#6A6A6A", linewidth=0.6))
    peaks = [0.20, 0.30, 0.42, 0.56, 0.68, 0.78]
    heights = [0.28, 0.46, 0.35, 0.62, 0.42, 0.30]
    for p, ht in zip(peaks, heights):
        ax.add_line(Line2D([x + w * p, x + w * p], [base_y, y + h * ht], color="#048BA8", linewidth=1.1))
    ax.add_patch(patches.Circle((x + w * 0.78, y + h * 0.68), w * 0.055, facecolor="#16DB93", edgecolor="#4B4B4B", linewidth=0.45))


def draw_gem_icon(ax, x: float, y: float, w: float, h: float) -> None:
    nodes = [
        (0.25, 0.55),
        (0.42, 0.72),
        (0.58, 0.48),
        (0.73, 0.66),
        (0.43, 0.32),
        (0.67, 0.28),
    ]
    edges = [(0, 1), (1, 2), (2, 3), (2, 4), (4, 5), (3, 5), (0, 4)]
    for a, b in edges:
        ax.add_line(
            Line2D(
                [x + w * nodes[a][0], x + w * nodes[b][0]],
                [y + h * nodes[a][1], y + h * nodes[b][1]],
                color="#777777",
                linewidth=0.7,
            )
        )
    for i, (nx, ny) in enumerate(nodes):
        color = "#4B7FAD" if i % 2 else "#9E4A4F"
        ax.add_patch(patches.Circle((x + w * nx, y + h * ny), w * 0.045, facecolor=color, edgecolor="white", linewidth=0.5))


def draw_heatmap_icon(ax, x: float, y: float, w: float, h: float) -> None:
    colors = ["#FFFFFF", "#D9ECEA", "#6FB7AE", "#104547"]
    for row in range(4):
        for col in range(5):
            idx = (row * 2 + col) % len(colors)
            ax.add_patch(
                patches.Rectangle(
                    (x + w * (0.15 + col * 0.13), y + h * (0.2 + row * 0.13)),
                    w * 0.12,
                    h * 0.12,
                    facecolor=colors[idx],
                    edgecolor="#EEEEEE",
                    linewidth=0.25,
                )
            )


def draw_exchange_icon(ax, x: float, y: float, w: float, h: float) -> None:
    ax.add_patch(patches.Circle((x + w * 0.30, y + h * 0.50), w * 0.10, facecolor=PHYLUM_COLORS["Actinobacteria"], edgecolor="white", linewidth=0.6))
    ax.add_patch(patches.Circle((x + w * 0.72, y + h * 0.50), w * 0.10, facecolor=PHYLUM_COLORS["Gammaproteobacteria"], edgecolor="white", linewidth=0.6))
    arrow(ax, (x + w * 0.39, y + h * 0.54), (x + w * 0.62, y + h * 0.54), color="#4B7FAD", lw=0.8, mutation_scale=7)
    arrow(ax, (x + w * 0.63, y + h * 0.44), (x + w * 0.40, y + h * 0.44), color="#9E4A4F", lw=0.8, mutation_scale=7)


def draw_network_icon(ax, x: float, y: float, w: float, h: float) -> None:
    pts = [(0.25, 0.38), (0.42, 0.66), (0.62, 0.44), (0.75, 0.70), (0.70, 0.25)]
    for a, b in [(0, 1), (1, 2), (2, 3), (2, 4), (0, 4)]:
        ax.add_line(Line2D([x + w * pts[a][0], x + w * pts[b][0]], [y + h * pts[a][1], y + h * pts[b][1]], color="#777777", linewidth=0.7))
    for idx, (px, py) in enumerate(pts):
        ax.add_patch(patches.Circle((x + w * px, y + h * py), w * 0.045, facecolor=list(PHYLUM_COLORS.values())[idx], edgecolor="white", linewidth=0.5))


def small_input(ax, x: float, y: float, w: float, h: float, title: str, icon: str) -> None:
    rounded_box(ax, x, y, w, h, facecolor=INPUT_FILL, edgecolor="#BFC8CC", linewidth=0.6, radius=0.45)
    icon_x, icon_y, icon_w, icon_h = x + 0.7, y + 0.8, w * 0.33, h - 1.6
    if icon == "genome":
        draw_genome_icon(ax, icon_x, icon_y, icon_w, icon_h)
    elif icon == "plate":
        draw_plate_icon(ax, icon_x, icon_y, icon_w, icon_h)
    elif icon == "metabolomics":
        draw_metabolomics_icon(ax, icon_x, icon_y, icon_w, icon_h)
    add_text(ax, x + w * 0.62, y + h * 0.55, title, size=7.2, weight="bold")


def analysis_card(ax, x: float, y: float, w: float, h: float, title: str, body: str, icon: str) -> None:
    rounded_box(ax, x, y, w, h, facecolor=ANALYSIS_FILL, edgecolor="#B8C8BD", linewidth=0.6, radius=0.5)
    icon_x, icon_y, icon_w, icon_h = x + 0.4, y + 0.7, w * 0.24, h - 1.4
    if icon == "heatmap":
        draw_heatmap_icon(ax, icon_x, icon_y, icon_w, icon_h)
    elif icon == "exchange":
        draw_exchange_icon(ax, icon_x, icon_y, icon_w, icon_h)
    elif icon == "network":
        draw_network_icon(ax, icon_x, icon_y, icon_w, icon_h)
    add_text(ax, x + w * 0.65, y + h * 0.66, title, size=7.4, weight="bold")
    add_text(ax, x + w * 0.65, y + h * 0.36, body, size=6.5, color=TEXT)


def output_card(ax, x: float, y: float, w: float, h: float, title: str, body: str) -> None:
    rounded_box(ax, x, y, w, h, facecolor=OUTPUT_FILL, edgecolor="#C7BACB", linewidth=0.6, radius=0.5)
    add_text(ax, x + w / 2, y + h * 0.66, title, size=7.2, weight="bold")
    add_text(ax, x + w / 2, y + h * 0.34, body, size=6.7, color=TEXT)


def main() -> None:
    apply_publication_style(
        extra_rc={
            "font.family": "sans-serif",
            "font.sans-serif": ["Helvetica", "Arial", "Avenir Next", "DejaVu Sans"],
            "font.size": 7,
            "axes.labelsize": 7,
            "axes.titlesize": 8,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "legend.fontsize": 7,
            "legend.title_fontsize": 7,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "savefig.bbox": None,
        }
    )

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=mm_to_in(180, 150))
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.axis("off")
    fig.subplots_adjust(left=0, right=1, top=1, bottom=0)

    add_text(
        ax,
        50,
        96.4,
        "Workflow linking SynCom traits, conditioned chemistry and model-predicted interactions",
        size=9.4,
        weight="bold",
    )

    # Section labels.
    for x, label in [
        (12.0, "Defined community"),
        (31.0, "Measured inputs"),
        (52.0, "Curated models"),
        (70.0, "Model analyses"),
        (90.0, "Ecological readouts"),
    ]:
        add_text(ax, x, 91.6, label, size=7.3, weight="bold", color=MUTED)
        ax.add_line(Line2D([x - 7.2, x + 7.2], [89.9, 89.9], color="#D0D0D0", linewidth=0.6))

    # Defined community.
    rounded_box(ax, 4, 62.5, 17.5, 24.0, facecolor="#FFFFFF", edgecolor="#BABABA", linewidth=0.65, radius=0.7)
    draw_syncom_icon(ax, 5.2, 68.0, 15.0, 12.0)
    add_text(ax, 12.75, 84.0, "Soil SynCom", size=8.0, weight="bold")
    add_text(ax, 12.75, 64.8, "taxonomically diverse\nisolate collection", size=7.0)

    # Input measurements.
    small_input(ax, 25.0, 76.2, 16.7, 9.0, "Genome\nsequences", "genome")
    small_input(ax, 25.0, 63.6, 16.7, 9.0, "Growth\nphenotypes", "plate")
    small_input(ax, 25.0, 51.0, 16.7, 9.0, "Basal and\nconditioned\nmetabolomics", "metabolomics")

    # Model reconstruction core.
    rounded_box(ax, 46.1, 62.0, 15.8, 24.5, facecolor=MODEL_FILL, edgecolor="#C7BBA9", linewidth=0.65, radius=0.7)
    draw_gem_icon(ax, 48.0, 72.5, 12.0, 8.2)
    add_text(ax, 54.0, 84.0, "Thermodynamic\nGEMs", size=8.0, weight="bold")
    add_text(ax, 54.0, 67.6, "draft reconstruction\nmanual curation\ngrowth validation", size=7.0)
    mini_arrow(ax, 48.2, 64.8, 59.8, 64.8)
    add_text(ax, 54.0, 61.1, "21 strain-resolved models", size=6.9, color=MUTED)

    # Environmental comparison strip feeding DiMEs.
    rounded_box(ax, 46.1, 47.8, 15.8, 8.4, facecolor=ENV_FILL, edgecolor="#B7C5D4", linewidth=0.6, radius=0.55)
    ax.add_patch(patches.Rectangle((48.1, 51.3), 5.1, 1.7, facecolor="#EFE6DC", edgecolor="#B7A999", linewidth=0.4))
    ax.add_patch(patches.Rectangle((54.8, 51.3), 5.1, 1.7, facecolor="#C9E3EC", edgecolor="#8EAFC0", linewidth=0.4))
    add_text(ax, 50.65, 49.7, "basal", size=6.7)
    add_text(ax, 57.35, 49.7, "conditioned", size=6.7)
    add_text(ax, 54.0, 54.3, "chemical contexts", size=7.2, weight="bold")

    # Analysis cards.
    analysis_card(
        ax,
        66.1,
        73.0,
        17.7,
        12.5,
        "iMM",
        "minimal\nnutrient sets",
        "heatmap",
    )
    analysis_card(
        ax,
        66.1,
        56.4,
        17.7,
        12.5,
        "DiMEs",
        "feasible uptake\nand secretion",
        "exchange",
    )
    analysis_card(
        ax,
        66.1,
        39.8,
        17.7,
        12.5,
        "NIS and\nessentiality",
        "biosynthetic\nroutes",
        "network",
    )

    # Output cards.
    output_card(ax, 86.7, 76.5, 11.0, 8.6, "Nutrient\ndependence", "core/variable\niMMs")
    output_card(ax, 86.7, 65.9, 11.0, 8.6, "Metabolic\nopportunity", "conditioned\nresource access")
    output_card(ax, 86.7, 55.3, 11.0, 8.6, "Niche\noverlap", "NOI\ncomparison")
    output_card(ax, 86.7, 44.7, 11.0, 8.6, "Candidate\nexchange", "producer-\nconsumer links")
    output_card(ax, 86.7, 34.1, 11.0, 8.6, "Biosynthetic\nroutes", "similarity and\nbottlenecks")

    rounded_box(ax, 8.0, 14.7, 84.5, 9.8, facecolor="#FFFFFF", edgecolor="#CFCFCF", linewidth=0.65, radius=0.65)
    add_text(
        ax,
        50.25,
        21.4,
        "Trait–environment interpretation",
        size=8.0,
        weight="bold",
    )
    add_text(
        ax,
        50.25,
        17.7,
        "Measured metabolomics defines the chemical context; curated GEMs identify model-predicted feasible\nnutrient requirements, exchange opportunities and candidate interactions.",
        size=7.0,
    )

    # Main arrows.
    arrow(ax, (21.5, 74.8), (24.8, 80.8), lw=0.8)
    arrow(ax, (21.5, 74.8), (24.8, 68.1), lw=0.8)
    arrow(ax, (21.5, 74.8), (24.8, 55.5), lw=0.8)
    arrow(ax, (41.9, 80.8), (45.9, 78.6), lw=0.8)
    arrow(ax, (41.9, 68.1), (45.9, 70.8), lw=0.8)
    arrow(ax, (41.9, 55.5), (45.9, 52.0), lw=0.8)
    arrow(ax, (61.9, 74.4), (65.9, 79.2), lw=0.8)
    arrow(ax, (61.9, 74.4), (65.9, 62.7), lw=0.8)
    arrow(ax, (61.9, 74.4), (65.9, 46.1), lw=0.8)
    arrow(ax, (61.9, 52.0), (65.9, 62.7), lw=0.8, color="#668090")

    arrow(ax, (83.8, 79.2), (86.5, 80.8), lw=0.8)
    arrow(ax, (83.8, 62.7), (86.5, 70.2), lw=0.8)
    arrow(ax, (83.8, 62.7), (86.5, 59.6), lw=0.8)
    arrow(ax, (83.8, 62.7), (86.5, 49.0), lw=0.8)
    arrow(ax, (83.8, 46.1), (86.5, 38.4), lw=0.8)

    # Phylum legend.
    legend_x, legend_y = 5.5, 5.8
    add_text(ax, legend_x, legend_y + 4.0, "Phylum colours", size=7.1, weight="bold", ha="left")
    lx = legend_x + 12.5
    for i, (name, color) in enumerate(PHYLUM_COLORS.items()):
        xx = lx + (i % 3) * 24.5
        yy = legend_y + 4.0 - (i // 3) * 3.5
        ax.add_patch(patches.Rectangle((xx, yy - 0.95), 1.6, 1.6, facecolor=color, edgecolor="none"))
        add_text(ax, xx + 2.2, yy - 0.15, name, size=6.9, ha="left")

    for ext in ("pdf", "svg", "png"):
        fig.savefig(f"{OUT_STEM}.{ext}", dpi=600, bbox_inches=fig.bbox_inches, pad_inches=0)
    plt.close(fig)


if __name__ == "__main__":
    main()
