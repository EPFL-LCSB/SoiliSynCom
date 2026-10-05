# File: src/pytfa_next/visualization/publication_style.py
"""Publication-quality matplotlib styling for pyTFA-next outputs.

Purpose
-------
One import, one call, no Illustrator. Use this module to make figures produced
by notebooks under ``scripts/`` or ``notebooks/`` land on the page already in
camera-ready form: consistent fonts, clean spines, outward ticks, legends
outside the plot area, and a muted paper-friendly palette.

The style targets Nature / Cell / PLOS-family conventions but is generic
enough for any high-impact journal. See ``docs/figure_style_guide.md`` for the
full specification.

Usage
-----

    from pytfa_next.visualization.publication_style import (
        apply_publication_style,
        setup_axes,
        add_panel_letter,
        save_publication_figure,
        PUBLICATION_COLORS,
    )

    apply_publication_style()  # call once at notebook start

    fig, ax = plt.subplots(figsize=(3.5, 3.5))
    ax.scatter(x, y, color=PUBLICATION_COLORS["substrate"])
    ax.set_xlabel("MOPhUs yield γᵢ = |αᵢ|")
    ax.set_ylabel("dTFA mean yield γᵢ = |vᵢ|/μ")
    setup_axes(ax)
    add_panel_letter(ax, "d")
    save_publication_figure(fig, "figures/panel_d", formats=("pdf", "png"))

Rationale
---------
- Helvetica → Arial fallback because most journals require Helvetica/Arial.
- ``pdf.fonttype = 42`` (TrueType) so the PDF remains editable in Illustrator
  as an escape hatch, and so text remains searchable.
- ``svg.fonttype = "none"`` so SVG output preserves text as text (not paths).
- Outward ticks and hidden top/right spines match modern Nature style.
"""
from __future__ import annotations

import logging
import warnings
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import matplotlib as mpl
import matplotlib.pyplot as plt

__all__ = [
    "PUBLICATION_COLORS",
    "FIGURE_SIZE_SINGLE",
    "FIGURE_SIZE_DOUBLE",
    "FIGURE_SIZE_PANEL",
    "FIGURE_SIZE_SQUARE_LARGE",
    "FONT_SIZE_AXIS_LABEL",
    "FONT_SIZE_TICK_LABEL",
    "FONT_SIZE_LEGEND",
    "FONT_SIZE_ANNOTATION",
    "FONT_SIZE_PANEL_LETTER",
    "apply_publication_style",
    "setup_axes",
    "add_panel_letter",
    "place_legend_outside",
    "color_secondary_axis",
    "add_identity_line_with_band",
    "label_points",
    "save_publication_figure",
    "publication_style_context",
]

logger = logging.getLogger(__name__)

# ============================================================================
# Palette — muted, print-friendly, legible in grayscale conversion.
# ============================================================================
# Convention used throughout pyTFA-next manuscripts:
#   substrate / uptake      → muted brick red
#   product / secretion     → muted blue
#   biomass / reference     → near-black (softer than #000)
#   accent / growth rate    → muted orange (dual-axis)
#   neutral / grouped minor → medium gray
#   reference line (y=x)    → dark gray
#   confidence band         → light gray (use alpha≈0.3)
# Hex values chosen to match the MOPhUs v5 figure palette so this style
# drops in without visual drift.

PUBLICATION_COLORS: dict[str, str] = {
    "substrate": "#9E4A4F",
    "product": "#4B7FAD",
    "biomass": "#1F1F1F",
    "accent_orange": "#E89B44",
    "neutral": "#9E9E9E",
    "reference_line": "#606060",
    "confidence_band": "#CCCCCC",
    # Extended palette for categorical data beyond the primary three
    "secondary_green": "#5D8A66",
    "secondary_purple": "#7A5A8C",
    "secondary_teal": "#3F8B8E",
}

# ============================================================================
# Standard figure dimensions (inches)
# Column widths: Nature single ≈ 89 mm; Nature double ≈ 183 mm.
# ============================================================================

FIGURE_SIZE_SINGLE: tuple[float, float] = (3.5, 3.0)      # ~89 mm × 76 mm
FIGURE_SIZE_DOUBLE: tuple[float, float] = (7.2, 3.5)      # ~183 mm × 89 mm
FIGURE_SIZE_PANEL: tuple[float, float] = (3.5, 3.5)       # square panel
FIGURE_SIZE_SQUARE_LARGE: tuple[float, float] = (5.0, 5.0)


# ============================================================================
# Font sizes — hard floor at 7 pt for axis labels per user mandate.
# ============================================================================

FONT_SIZE_AXIS_LABEL: int = 8
FONT_SIZE_TICK_LABEL: int = 7
FONT_SIZE_LEGEND: int = 7
FONT_SIZE_ANNOTATION: int = 6
FONT_SIZE_PANEL_LETTER: int = 12


# ============================================================================
# Line / tick / spine geometry
# ============================================================================

SPINE_WIDTH: float = 0.75
TICK_WIDTH: float = 0.75
TICK_LENGTH: float = 3.0
MINOR_TICK_LENGTH: float = 1.75
DATA_LINE_WIDTH: float = 1.0
DEFAULT_MARKER_SIZE: float = 5.0

# ============================================================================
# Main entry point
# ============================================================================


def apply_publication_style(
    *,
    font_family: str | Sequence[str] | None = None,
    extra_rc: Mapping[str, Any] | None = None,
) -> None:
    """Apply global matplotlib rcParams for publication-quality figures.

    Call once at notebook start (ideally immediately after ``import matplotlib``).
    Idempotent — safe to call multiple times.

    Parameters
    ----------
    font_family
        Override the sans-serif font chain. Accepts either a single name
        (``"Helvetica"``) or a list of fallbacks. Defaults to a chain that
        works on macOS, Linux and Windows without manual font installs.
    extra_rc
        Extra rcParams to merge on top (e.g. ``{"figure.dpi": 200}``).

    Notes
    -----
    If Helvetica is unavailable, matplotlib silently falls back to the next
    entry in the list. You can verify what's actually being used with::

        import matplotlib.font_manager as fm
        print(fm.findfont("Helvetica"))

    On Linux, install ``fonts-liberation`` (apt) or ``liberation-fonts-common``
    (dnf) — Liberation Sans is a metric-compatible Helvetica substitute.
    """
    if font_family is None:
        font_chain: list[str] = [
            "Helvetica",
            "Helvetica Neue",
            "Arial",
            "Liberation Sans",
            "DejaVu Sans",
        ]
    elif isinstance(font_family, str):
        font_chain = [font_family, "Arial", "Liberation Sans", "DejaVu Sans"]
    else:
        font_chain = list(font_family)

    rc: dict[str, Any] = {
        # ----- typography -----
        "font.family": "sans-serif",
        "font.sans-serif": font_chain,
        "font.size": FONT_SIZE_TICK_LABEL,
        "axes.labelsize": FONT_SIZE_AXIS_LABEL,
        "axes.titlesize": FONT_SIZE_AXIS_LABEL,
        "axes.labelweight": "normal",
        "axes.titleweight": "normal",
        "xtick.labelsize": FONT_SIZE_TICK_LABEL,
        "ytick.labelsize": FONT_SIZE_TICK_LABEL,
        "legend.fontsize": FONT_SIZE_LEGEND,
        "legend.title_fontsize": FONT_SIZE_LEGEND,
        # ----- axes / spines -----
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.spines.bottom": True,
        "axes.spines.left": True,
        "axes.linewidth": SPINE_WIDTH,
        "axes.edgecolor": "#1F1F1F",
        "axes.labelcolor": "#1F1F1F",
        "axes.axisbelow": True,  # gridlines behind data if ever enabled
        # ----- ticks -----
        "xtick.direction": "out",
        "ytick.direction": "out",
        "xtick.major.size": TICK_LENGTH,
        "ytick.major.size": TICK_LENGTH,
        "xtick.minor.size": MINOR_TICK_LENGTH,
        "ytick.minor.size": MINOR_TICK_LENGTH,
        "xtick.major.width": TICK_WIDTH,
        "ytick.major.width": TICK_WIDTH,
        "xtick.minor.width": TICK_WIDTH * 0.75,
        "ytick.minor.width": TICK_WIDTH * 0.75,
        "xtick.color": "#1F1F1F",
        "ytick.color": "#1F1F1F",
        "xtick.major.pad": 3.0,
        "ytick.major.pad": 3.0,
        # ----- lines / markers -----
        "lines.linewidth": DATA_LINE_WIDTH,
        "lines.markersize": DEFAULT_MARKER_SIZE,
        "lines.markeredgewidth": 0.6,
        # ----- error bars (Nature style: capless, thin) -----
        "errorbar.capsize": 0.0,
        # ----- grid (off by default) -----
        "axes.grid": False,
        "grid.linewidth": 0.4,
        "grid.color": "#E5E5E5",
        # ----- legend -----
        "legend.frameon": False,
        "legend.handlelength": 1.5,
        "legend.handletextpad": 0.5,
        "legend.borderaxespad": 0.3,
        "legend.columnspacing": 1.0,
        # ----- figure / saving -----
        "figure.dpi": 150,
        "savefig.dpi": 600,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.05,
        "savefig.transparent": False,
        "pdf.fonttype": 42,   # TrueType: editable in Illustrator, searchable text
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        # ----- math text: let mathtext use the same sans-serif family -----
        "mathtext.fontset": "dejavusans",
        "mathtext.default": "regular",
    }

    if extra_rc:
        rc.update(extra_rc)

    mpl.rcParams.update(rc)
    logger.info("Applied pyTFA-next publication figure style (font chain=%s).", font_chain)


# ============================================================================
# Per-axis helpers
# ============================================================================


def setup_axes(
    ax,
    *,
    hide_top: bool = True,
    hide_right: bool = True,
    hide_left: bool = False,
    hide_bottom: bool = False,
) -> None:
    """Apply per-axis publication styling.

    Usually redundant after :func:`apply_publication_style`, but needed when:
    - an axis was created before the style was applied,
    - a ``twinx`` axis needs the right spine re-shown,
    - you want a more minimal axis (e.g. hide both left and bottom for inset).
    """
    ax.spines["top"].set_visible(not hide_top)
    ax.spines["right"].set_visible(not hide_right)
    ax.spines["left"].set_visible(not hide_left)
    ax.spines["bottom"].set_visible(not hide_bottom)
    for spine in ax.spines.values():
        spine.set_linewidth(SPINE_WIDTH)
        spine.set_color("#1F1F1F")
    ax.tick_params(
        axis="both",
        direction="out",
        width=TICK_WIDTH,
        length=TICK_LENGTH,
        color="#1F1F1F",
    )


def add_panel_letter(
    ax,
    letter: str,
    *,
    x_offset: float = -0.15,
    y_offset: float = 1.02,
    fontsize: int = FONT_SIZE_PANEL_LETTER,
    fontweight: str = "bold",
) -> Any:
    """Add a bold panel letter (a, b, c, ...) in the top-left, Nature-style.

    The default position sits just above the top of the axes and slightly
    to the left of the y-axis label. For multi-row figures you usually want
    to tune ``x_offset`` per panel so letters visually align.
    """
    return ax.text(
        x_offset,
        y_offset,
        letter,
        transform=ax.transAxes,
        fontsize=fontsize,
        fontweight=fontweight,
        va="bottom",
        ha="left",
    )


def place_legend_outside(
    ax,
    *,
    loc: str = "right",
    pad: float = 0.02,
    ncol: int = 1,
    **legend_kwargs: Any,
) -> Any:
    """Place the legend outside the plot area.

    Parameters
    ----------
    loc : {'right', 'top', 'bottom', 'top-right'}
        Which side (or corner) of the axes to anchor the legend to.
    pad : float
        Fractional offset from the axes edge (axes coordinates).
    ncol : int
        Number of legend columns (useful for ``loc='top'`` or ``'bottom'``).
    **legend_kwargs
        Forwarded to ``ax.legend``.
    """
    configs = {
        "right": dict(bbox_to_anchor=(1.0 + pad, 1.0), loc="upper left"),
        "top-right": dict(bbox_to_anchor=(1.0 + pad, 1.0), loc="upper left"),
        "top": dict(bbox_to_anchor=(0.0, 1.0 + pad), loc="lower left"),
        "top-center": dict(bbox_to_anchor=(0.5, 1.0 + pad), loc="lower center"),
        "bottom": dict(bbox_to_anchor=(0.0, -0.15 - pad), loc="upper left"),
        "bottom-center": dict(bbox_to_anchor=(0.5, -0.15 - pad), loc="upper center"),
    }
    if loc not in configs:
        raise ValueError(
            f"loc must be one of {list(configs)}, got {loc!r}"
        )
    cfg = configs[loc]
    return ax.legend(
        bbox_to_anchor=cfg["bbox_to_anchor"],
        loc=cfg["loc"],
        frameon=False,
        ncol=ncol,
        **legend_kwargs,
    )


def color_secondary_axis(ax_twin, color: str) -> None:
    """Color a ``twinx`` secondary axis (spine, ticks, label) uniformly.

    Matches the Nature dual-axis convention used in Fig. 4a of the MOPhUs paper
    (orange right axis for growth rate on top of the black left axis for yield).
    """
    ax_twin.spines["right"].set_visible(True)
    ax_twin.spines["right"].set_color(color)
    ax_twin.spines["right"].set_linewidth(SPINE_WIDTH)
    ax_twin.tick_params(axis="y", colors=color, width=TICK_WIDTH, length=TICK_LENGTH)
    ax_twin.yaxis.label.set_color(color)


def add_identity_line_with_band(
    ax,
    *,
    band_fraction: float = 0.10,
    line_color: str = PUBLICATION_COLORS["reference_line"],
    band_color: str = PUBLICATION_COLORS["confidence_band"],
    line_kwargs: Mapping[str, Any] | None = None,
    band_kwargs: Mapping[str, Any] | None = None,
    label: str = "y = x",
    band_label: str | None = None,
) -> tuple[Any, Any]:
    """Add a dashed identity line and a ±fraction shaded band (Fig 4d style).

    Operates in axes (display) coordinates by reading the current xlim/ylim,
    so call after plotting the data.

    Parameters
    ----------
    band_fraction : float
        Half-width of the tolerance band as a fraction (``0.10`` → ±10%).
        Set to 0 to disable the band.
    label, band_label : str
        Legend labels; ``band_label`` defaults to ``f"±{band_fraction*100:.0f}%"``.
    """
    import numpy as np

    if band_label is None:
        band_label = f"±{band_fraction * 100:.0f}%" if band_fraction > 0 else None

    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    lo = min(x0, y0)
    hi = max(x1, y1)
    xs = np.geomspace(max(lo, 1e-30), hi, 200) if lo > 0 and ax.get_xscale() == "log" else np.linspace(lo, hi, 200)

    line_defaults = dict(linestyle="--", linewidth=0.75, zorder=1)
    line_defaults.update(line_kwargs or {})
    (line_artist,) = ax.plot(xs, xs, color=line_color, label=label, **line_defaults)

    band_artist = None
    if band_fraction > 0:
        band_defaults = dict(alpha=0.3, linewidth=0, zorder=0)
        band_defaults.update(band_kwargs or {})
        band_artist = ax.fill_between(
            xs,
            xs * (1.0 - band_fraction),
            xs * (1.0 + band_fraction),
            color=band_color,
            label=band_label,
            **band_defaults,
        )

    # Restore the original data limits (plot() above may have stretched them).
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    return line_artist, band_artist


def label_points(
    ax,
    xs: Iterable[float],
    ys: Iterable[float],
    labels: Iterable[str],
    *,
    offset: tuple[float, float] = (4.0, 4.0),
    fontsize: int = FONT_SIZE_ANNOTATION,
    color: str = PUBLICATION_COLORS["neutral"],
    **text_kwargs: Any,
) -> list[Any]:
    """Annotate each point with a small label offset in display points."""
    artists = []
    for x, y, lab in zip(xs, ys, labels):
        art = ax.annotate(
            lab,
            (x, y),
            xytext=offset,
            textcoords="offset points",
            fontsize=fontsize,
            color=color,
            **text_kwargs,
        )
        artists.append(art)
    return artists


# ============================================================================
# Saving
# ============================================================================


def save_publication_figure(
    fig,
    output_stem: Path | str,
    *,
    formats: Sequence[str] = ("pdf", "png"),
    dpi: int = 600,
    overwrite: bool = True,
    transparent: bool = False,
) -> list[Path]:
    """Save a figure in multiple formats with publication-quality settings.

    Parameters
    ----------
    output_stem
        File path WITHOUT extension. Each format in ``formats`` gets appended.
    formats
        File extensions, e.g. ``("pdf", "png")`` or ``("pdf", "svg", "png")``.
    dpi
        Raster DPI (applies to PNG/JPEG/TIFF; PDF/SVG are vector).
    transparent
        True for transparent background (useful for overlay on slides).

    Returns
    -------
    list of paths written.
    """
    stem = Path(output_stem)
    stem.parent.mkdir(parents=True, exist_ok=True)

    written: list[Path] = []
    for fmt in formats:
        path = stem.with_suffix(f".{fmt}")
        if path.exists() and not overwrite:
            warnings.warn(f"Skipping existing figure: {path}", stacklevel=2)
            continue
        fig.savefig(
            path,
            dpi=dpi,
            bbox_inches="tight",
            pad_inches=0.05,
            transparent=transparent,
        )
        written.append(path)
        logger.info("Wrote %s", path)
    return written


# ============================================================================
# Context manager
# ============================================================================


@contextmanager
def publication_style_context(overrides: Mapping[str, Any] | None = None):
    """Temporarily apply the publication style within a ``with`` block.

    Named ``publication_style_context`` (not ``publication_style``) to avoid a
    name collision with the ``publication_style`` submodule attribute on the
    parent ``pytfa_next.visualization`` package — re-exporting a function under
    the same name as its parent module shadows the module on the package and
    breaks ``import pytfa_next.visualization.publication_style as ps``.

    Usage::

        with publication_style_context():
            fig, ax = plt.subplots()
            save_publication_figure(fig, "my_figure")

        # with custom overrides (use rcParam-style dotted keys):
        with publication_style_context({"figure.dpi": 200}):
            ...

    The caller's previous rcParams are restored on exit. ``overrides`` must be
    a dict because rcParam keys (e.g. ``"figure.dpi"``) contain dots and can't
    be expressed as Python keyword arguments.
    """
    original = dict(mpl.rcParams)
    try:
        apply_publication_style(extra_rc=overrides)
        yield
    finally:
        mpl.rcParams.update(original)
