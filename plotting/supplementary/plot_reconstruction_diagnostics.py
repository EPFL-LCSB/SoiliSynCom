from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
PROJECT_ROOT = HERE.parents[1]
SOURCE = PROJECT_ROOT / "tables" / "graphing" / "GEMStats" / "GEMStats.xlsx"
OUTPUT_STEM = PROJECT_ROOT / "outputs" / "reconstruction_diagnostics" / "S05_diagnostics"

sys.path.insert(0, str(PROJECT_ROOT))
from publication_style import (  # noqa: E402
    PUBLICATION_COLORS,
    add_panel_letter,
    apply_publication_style,
    setup_axes,
)


NATURE_2COL_WIDTH_IN = 180 / 25.4
NATURE_FIGURE_HEIGHT_IN = 74 / 25.4


GROUP_ORDER = [
    "Alphaproteobacteria",
    "Betaproteobacteria",
    "Gammaproteobacteria",
    "Actinobacteria",
    "Bacteroidetes",
    "Firmicutes",
]

GROUP_COLORS = {
    "Alphaproteobacteria": "#3B6EA8",
    "Betaproteobacteria": "#5BA7C9",
    "Gammaproteobacteria": PUBLICATION_COLORS["secondary_green"],
    "Actinobacteria": PUBLICATION_COLORS["substrate"],
    "Bacteroidetes": PUBLICATION_COLORS["secondary_purple"],
    "Firmicutes": PUBLICATION_COLORS["neutral"],
}


def linear_fit_stats(x: pd.Series, y: pd.Series) -> tuple[np.ndarray, np.ndarray, float, float, float]:
    x_arr = x.to_numpy(dtype=float)
    y_arr = y.to_numpy(dtype=float)
    slope, intercept = np.polyfit(x_arr, y_arr, 1)
    fit_x = np.linspace(x_arr.min() - 0.2, x_arr.max() + 0.2, 200)
    fit_y = slope * fit_x + intercept
    r = float(np.corrcoef(x_arr, y_arr)[0, 1])
    residual_sd = float(np.std(y_arr - (slope * x_arr + intercept)))
    return fit_x, fit_y, r, r**2, residual_sd


def add_taxonomic_group(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["Phylum"] = df["Phylum"].str.strip()
    df["Taxonomic group"] = df["Phylum"]
    proteo = df["Phylum"] == "Proteobacteria"
    df.loc[proteo, "Taxonomic group"] = df.loc[proteo, "Class"]
    return df


def scatter_by_group(ax, df: pd.DataFrame, x_col: str, y_col: str) -> list:
    handles = []
    for group in GROUP_ORDER:
        sub = df[df["Taxonomic group"] == group]
        if sub.empty:
            continue
        artist = ax.scatter(
            sub[x_col],
            sub[y_col],
            s=33,
            color=GROUP_COLORS[group],
            edgecolor="white",
            linewidth=0.55,
            label=group,
            zorder=3,
        )
        handles.append(artist)
    return handles


def annotate_fit(ax, x: pd.Series, y: pd.Series, loc: tuple[float, float] = (0.04, 0.95)) -> None:
    fit_x, fit_y, r, r2, residual_sd = linear_fit_stats(x, y)
    ax.plot(fit_x, fit_y, color=PUBLICATION_COLORS["reference_line"], linewidth=0.9, zorder=2)
    ax.fill_between(
        fit_x,
        fit_y - residual_sd,
        fit_y + residual_sd,
        color=PUBLICATION_COLORS["confidence_band"],
        alpha=0.32,
        linewidth=0,
        zorder=1,
    )
    ax.text(
        loc[0],
        loc[1],
        f"r = {r:.2f}\nR$^2$ = {r2:.2f}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=6,
        bbox={
            "boxstyle": "round,pad=0.22",
            "facecolor": "white",
            "edgecolor": "#D7D7D7",
            "linewidth": 0.45,
        },
    )


def plot_size_relationships(ax, df: pd.DataFrame, y_col: str, ylabel: str, title: str) -> None:
    scatter_by_group(ax, df, "Total Length (Mb)", y_col)
    annotate_fit(ax, df["Total Length (Mb)"], df[y_col])
    ax.set_xlabel("Genome size (Mb)")
    ax.set_ylabel(ylabel)
    ax.set_title(title, pad=4)
    ax.set_xlim(df["Total Length (Mb)"].min() - 0.45, df["Total Length (Mb)"].max() + 0.45)
    pad = (df[y_col].max() - df[y_col].min()) * 0.13
    ax.set_ylim(df[y_col].min() - pad, df[y_col].max() + pad)
    setup_axes(ax)


def plot_validation_improvement(ax, df: pd.DataFrame) -> None:
    metric_specs = [
        ("F1", "F1Draft", "F1Final", 0.0, PUBLICATION_COLORS["product"]),
        ("MCC", "MCCDraft", "MCCFinal", 2.7, PUBLICATION_COLORS["substrate"]),
    ]

    for label, draft_col, final_col, base_x, color in metric_specs:
        x_draft = base_x
        x_final = base_x + 0.85
        for _, row in df.iterrows():
            ax.plot(
                [x_draft, x_final],
                [row[draft_col], row[final_col]],
                color="#B7B7B7",
                linewidth=0.55,
                alpha=0.65,
                zorder=1,
            )
            ax.scatter(
                [x_draft, x_final],
                [row[draft_col], row[final_col]],
                s=10,
                color="#D0D0D0",
                edgecolor="white",
                linewidth=0.25,
                zorder=2,
            )

        means = [df[draft_col].mean(), df[final_col].mean()]
        ax.plot([x_draft, x_final], means, color=color, linewidth=1.45, zorder=4)
        ax.scatter(
            [x_draft, x_final],
            means,
            s=34,
            color=color,
            edgecolor="white",
            linewidth=0.7,
            zorder=5,
        )
        ax.text(
            (x_draft + x_final) / 2,
            1.06,
            label,
            ha="center",
            va="bottom",
            fontsize=7,
        )
        ax.text(
            x_draft,
            means[0] - 0.08,
            f"{means[0]:.2f}",
            ha="center",
            va="top",
            fontsize=6,
            color=color,
        )
        ax.text(
            x_final,
            min(means[1] + 0.06, 1.03),
            f"{means[1]:.2f}",
            ha="center",
            va="bottom",
            fontsize=6,
            color=color,
        )

    ax.axhline(0, color="#CFCFCF", linewidth=0.6, zorder=0)
    ax.set_xticks([0.0, 0.85, 2.7, 3.55])
    ax.set_xticklabels(["Draft", "Curated", "Draft", "Curated"])
    ax.set_ylabel("Validation score")
    ax.set_title("Curation improves model validation", pad=4)
    ax.set_xlim(-0.45, 4.0)
    ax.set_ylim(-0.38, 1.12)
    setup_axes(ax)


def save_nature_figure(fig, output_stem: Path, formats: tuple[str, ...] = ("pdf", "png", "svg")) -> list[Path]:
    output_stem.parent.mkdir(parents=True, exist_ok=True)
    written = []
    for fmt in formats:
        path = output_stem.with_suffix(f".{fmt}")
        with plt.rc_context({"savefig.bbox": "standard", "savefig.pad_inches": 0.0}):
            fig.savefig(path, dpi=600, bbox_inches=None, pad_inches=0, transparent=False)
        written.append(path)
    return written


def main() -> None:
    apply_publication_style(
        extra_rc={
            "font.size": 6,
            "axes.labelsize": 7,
            "axes.titlesize": 7,
            "xtick.labelsize": 6,
            "ytick.labelsize": 6,
            "legend.fontsize": 6,
            "legend.title_fontsize": 6,
            "savefig.bbox": "standard",
            "savefig.pad_inches": 0.0,
        }
    )
    stats = pd.read_excel(SOURCE, sheet_name="Sheet2")
    stats = add_taxonomic_group(stats)
    validation = pd.read_excel(SOURCE, sheet_name="Sheet4")

    fig, axes = plt.subplots(1, 3, figsize=(NATURE_2COL_WIDTH_IN, NATURE_FIGURE_HEIGHT_IN))
    ax_a, ax_b, ax_c = axes.ravel()

    plot_size_relationships(
        ax_a,
        stats,
        y_col="Rxns",
        ylabel="Curated model reactions",
        title="Genome size versus model reactions",
    )
    plot_size_relationships(
        ax_b,
        stats,
        y_col="Genes",
        ylabel="Model-associated metabolic genes",
        title="Genome size versus metabolic genes",
    )
    plot_validation_improvement(ax_c, validation)

    for letter, ax in zip("abc", axes.ravel(), strict=True):
        add_panel_letter(ax, letter, x_offset=-0.18, y_offset=1.04, fontsize=7)

    handles, labels = ax_a.get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        title="Taxonomic group",
        loc="lower center",
        bbox_to_anchor=(0.5, 0.02),
        ncol=3,
        frameon=False,
        handletextpad=0.4,
        columnspacing=0.9,
    )
    fig.subplots_adjust(left=0.075, right=0.99, top=0.84, bottom=0.36, wspace=0.42)

    written = save_nature_figure(fig, OUTPUT_STEM)
    for path in written:
        print(path)


if __name__ == "__main__":
    main()
