"""droplet_statistics.py — Stage 3: per-card statistics from detected droplets.

For each card this module computes:
  - Covered area (%)
  - Droplets per card, Droplets per cm²
  - Diameter statistics, both raw (stain diameter) and spread-factor-corrected
    (real drop diameter):
      NMD, NMD10, NMD90  — numeric percentile diameters
      VMD, VMD10, VMD90  — volume percentile diameters
      CH                  — Coefficient of Homogeneity  (VMD / NMD)
      Mean diameter, Diameter STD

Results are written to a CSV file and also returned as a pandas DataFrame.

Spread factor
-------------
The spread factor relates the stain diameter visible on the WSP paper to the
actual drop diameter that produced it.  The correction polynomial was fitted
from calibration data (see MATLAB spread_factor_equation.m):

    stain_diams  = [100, 200, 300, 400, 500, 600]  µm
    drop_diams   = [ 59, 109, 155, 200, 243, 285]  µm

A linear-in-log regression gives:
    drop_diam = a + b * ln(stain_diam)      (a, b fit from data)

which in the MATLAB code was approximated as:
    drop_diam = 0.53549306 * stain_diam - 0.000084839 * stain_diam²

Both forms are available; this code uses the MATLAB polynomial to preserve
numerical parity with existing results.

Configuration
-------------
Card area and DPI are read from ``config.yaml`` via the shared ``Config``
dataclass from card_masks.py.

Usage
-----
    python droplet_statistics.py <image_file> [--output <dir>] [--config config.yaml]
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from card_masks import (
    Config, load_config,
    _build_card_mask, _cleanup, _get_regions, _order_cards,
    load_image,
)
from detect_droplets import detect_droplets, Droplet


# Physical card dimensions (standard WSP card)
_CARD_WIDTH_CM  = 2.6
_CARD_HEIGHT_CM = 7.6
_CARD_AREA_CM2  = _CARD_WIDTH_CM * _CARD_HEIGHT_CM


# ── Spread factor correction ──────────────────────────────────────────────────

def apply_spread_factor(stain_diam_um: np.ndarray) -> np.ndarray:
    """Convert stain diameters (µm) to real drop diameters (µm).

    Polynomial fitted from WSP calibration data (matches MATLAB implementation):
        drop_diam = 0.53549306 * stain_diam - 0.000084839 * stain_diam²
    """
    d = stain_diam_um
    return 0.53549306 * d - 0.000084839 * d**2


# ── Per-card statistics ───────────────────────────────────────────────────────

@dataclass
class CardStats:
    card_name:        str
    wsp_id:           str

    # Coverage
    covered_area_pct: float
    droplet_count:    int
    droplets_per_cm2: float

    # Diameter statistics — raw stain diameter
    nmd:          float
    nmd10:        float
    nmd90:        float
    vmd:          float
    vmd10:        float
    vmd90:        float
    ch:           float
    mean_diam:    float
    std_diam:     float

    # Diameter statistics — spread-factor corrected
    nmd_sf:       float
    nmd10_sf:     float
    nmd90_sf:     float
    vmd_sf:       float
    vmd10_sf:     float
    vmd90_sf:     float
    ch_sf:        float
    mean_diam_sf: float
    std_diam_sf:  float


def compute_card_stats(
    card_name: str,
    wsp_id: str,
    droplets: list[Droplet],
    mask_area_px: int,
    stain_area_px: int,
    dpi: int,
) -> CardStats:
    """Compute all statistics for one card.

    Parameters
    ----------
    card_name:
        Scan file name (used as an identifier in the output table).
    wsp_id:
        Card label, e.g. "WSP-1".
    droplets:
        List of Droplet objects from Stage 2.
    mask_area_px:
        Total number of pixels in the card mask — used as coverage denominator.
        Matches MATLAB's ``sum(maschere_cartina{idx}(:))``.
    stain_area_px:
        Raw stain pixel count before sub-segmentation (``R < threshold`` pixels).
        Matches MATLAB's ``sum(bw(:))`` for covered-area computation.
    dpi:
        Scanner resolution; determines the physical pixel size.
    """
    # Guard: no droplets detected
    nan = float("nan")
    if not droplets:
        return CardStats(
            card_name=card_name, wsp_id=wsp_id,
            covered_area_pct=0.0, droplet_count=0, droplets_per_cm2=0.0,
            nmd=nan, nmd10=nan, nmd90=nan,
            vmd=nan, vmd10=nan, vmd90=nan,
            ch=nan, mean_diam=nan, std_diam=nan,
            nmd_sf=nan, nmd10_sf=nan, nmd90_sf=nan,
            vmd_sf=nan, vmd10_sf=nan, vmd90_sf=nan,
            ch_sf=nan, mean_diam_sf=nan, std_diam_sf=nan,
        )

    pixel_um = 25_400.0 / dpi          # microns per pixel at given DPI

    # Build per-pixel-area count histogram (area_hist[k] = number of droplets with area k px)
    areas     = np.array([d.area_px for d in droplets], dtype=int)
    max_area  = int(areas.max())
    area_hist = np.bincount(areas, minlength=max_area + 1)[1:]  # index 0 unused; shift by 1

    # Indices 1..max_area map to pixel counts; compute equivalent stain diameters
    idx = np.arange(1, len(area_hist) + 1, dtype=float)
    stain_diam = pixel_um * np.sqrt(4.0 * idx / np.pi)

    # Coverage — Bug 2 fix: use raw stain pixels (R<threshold), not sum of subseg areas.
    # Matches MATLAB: 100 * sum(bw(:)) / sum(mask(:))
    covered_area_pct = (
        100.0 * stain_area_px / mask_area_px if mask_area_px > 0 else 0.0
    )

    # Density — Bug 3 fix: always use the fixed physical card area (2.6×7.6 cm²).
    # Matches MATLAB statistiche_gocce2: area_cartina = 2.6*7.6
    droplet_count    = int(area_hist.sum())
    droplets_per_cm2 = droplet_count / _CARD_AREA_CM2

    # Compute diameter statistics for both raw and spread-factor-corrected diameters
    raw_stats = _diameter_stats(stain_diam, area_hist)
    sf_diam   = apply_spread_factor(stain_diam)
    sf_stats  = _diameter_stats(sf_diam, area_hist)

    return CardStats(
        card_name=card_name, wsp_id=wsp_id,
        covered_area_pct=covered_area_pct,
        droplet_count=droplet_count,
        droplets_per_cm2=droplets_per_cm2,
        nmd=raw_stats["nmd"],   nmd10=raw_stats["p10"],  nmd90=raw_stats["p90"],
        vmd=raw_stats["vmd"],   vmd10=raw_stats["v10"],  vmd90=raw_stats["v90"],
        ch=raw_stats["ch"],
        mean_diam=raw_stats["mean"], std_diam=raw_stats["std"],
        nmd_sf=sf_stats["nmd"],  nmd10_sf=sf_stats["p10"],  nmd90_sf=sf_stats["p90"],
        vmd_sf=sf_stats["vmd"],  vmd10_sf=sf_stats["v10"],  vmd90_sf=sf_stats["v90"],
        ch_sf=sf_stats["ch"],
        mean_diam_sf=sf_stats["mean"], std_diam_sf=sf_stats["std"],
    )


def _diameter_stats(diameters: np.ndarray, counts: np.ndarray) -> dict:
    """Compute NMD/VMD family and mean/std from a diameter array and count weights.

    Parameters
    ----------
    diameters:
        1-D array of diameter values (one per pixel-area class).
    counts:
        1-D array of droplet counts for each diameter class (same length).
    """
    # Expand to per-droplet list for mean/std (weighted)
    total = float(counts.sum())
    mean  = float(np.average(diameters, weights=counts))
    # Bug 4 fix: N-1 (sample) variance — matches MATLAB std() default
    var   = (
        float(np.average((diameters - mean)**2, weights=counts) * total / (total - 1))
        if total > 1 else 0.0
    )
    std   = float(np.sqrt(max(var, 0.0)))

    # NMD: percentile diameters by cumulative count
    cum_n = np.cumsum(counts)
    nmd   = _weighted_percentile(diameters, cum_n, total, 0.50)
    p10   = _weighted_percentile(diameters, cum_n, total, 0.10)
    p90   = _weighted_percentile(diameters, cum_n, total, 0.90)

    # VMD: percentile diameters by cumulative volume (∝ d³)
    volumes = counts * ((4.0 / 3.0) * np.pi * (diameters / 2.0)**3)
    cum_v   = np.cumsum(volumes)
    total_v = cum_v[-1]
    vmd  = _weighted_percentile(diameters, cum_v, total_v, 0.50)
    v10  = _weighted_percentile(diameters, cum_v, total_v, 0.10)
    v90  = _weighted_percentile(diameters, cum_v, total_v, 0.90)

    ch   = vmd / nmd if nmd > 0 else float("nan")

    return dict(nmd=nmd, p10=p10, p90=p90,
                vmd=vmd, v10=v10, v90=v90,
                ch=ch, mean=mean, std=std)


def _weighted_percentile(
    values: np.ndarray,
    cumulative: np.ndarray,
    total: float,
    fraction: float,
) -> float:
    """Return the first value in *values* where *cumulative* exceeds fraction*total.

    Matches the MATLAB index arithmetic:
        idx = 1 + sum(cum <= fraction * total)
    """
    threshold = fraction * total
    idx = int(np.sum(cumulative <= threshold))
    idx = min(idx, len(values) - 1)
    return float(values[idx])


# ── Output ────────────────────────────────────────────────────────────────────

_CSV_COLUMNS = [
    "card_name", "WSP_id",
    "NMD", "NMD10", "NMD90", "VMD", "VMD10", "VMD90", "CH",
    "Mean diameter", "Diameter STD",
    "NMD_SF", "NMD10_SF", "NMD90_SF", "VMD_SF", "VMD10_SF", "VMD90_SF", "CH_SF",
    "Mean diameter_SF", "Diameter STD_SF",
    "Covered area", "Droplets per card", "Droplets per cm2",
]


def stats_to_dataframe(stats_list: list[CardStats]) -> pd.DataFrame:
    rows = []
    for s in stats_list:
        rows.append({
            "card_name":         s.card_name,
            "WSP_id":            s.wsp_id,
            "NMD":               s.nmd,
            "NMD10":             s.nmd10,
            "NMD90":             s.nmd90,
            "VMD":               s.vmd,
            "VMD10":             s.vmd10,
            "VMD90":             s.vmd90,
            "CH":                s.ch,
            "Mean diameter":     s.mean_diam,
            "Diameter STD":      s.std_diam,
            "NMD_SF":            s.nmd_sf,
            "NMD10_SF":          s.nmd10_sf,
            "NMD90_SF":          s.nmd90_sf,
            "VMD_SF":            s.vmd_sf,
            "VMD10_SF":          s.vmd10_sf,
            "VMD90_SF":          s.vmd90_sf,
            "CH_SF":             s.ch_sf,
            "Mean diameter_SF":  s.mean_diam_sf,
            "Diameter STD_SF":   s.std_diam_sf,
            "Covered area":      s.covered_area_pct,
            "Droplets per card": s.droplet_count,
            "Droplets per cm2":  s.droplets_per_cm2,
        })
    return pd.DataFrame(rows, columns=_CSV_COLUMNS)


def print_stats_table(stats_list: list[CardStats]) -> None:
    """Print a formatted statistics table to stdout."""
    if not stats_list:
        return

    lbl_w = 19   # label column width (fits "Droplets per card" = 17 chars)
    val_w = 9    # value column width

    def _hline():
        print("-" * (lbl_w + val_w * len(stats_list)))

    def _hdr(label, values):
        print(f"{label:<{lbl_w}}" + "".join(f"{v:>{val_w}}" for v in values))

    def _num(label, values):
        cells = []
        for v in values:
            cells.append(f"{'N/A':>{val_w}}" if (isinstance(v, float) and np.isnan(v))
                         else f"{v:{val_w}.2f}")
        print(f"{label:<{lbl_w}}" + "".join(cells))

    def _pct(label, values):
        cells = []
        for v in values:
            cells.append(f"{'N/A':>{val_w}}" if (isinstance(v, float) and np.isnan(v))
                         else f"{v:{val_w-1}.2f}%")
        print(f"{label:<{lbl_w}}" + "".join(cells))

    def _int(label, values):
        print(f"{label:<{lbl_w}}" + "".join(f"{int(v):{val_w}d}" for v in values))

    print()
    _hdr("Card analysis:", [s.wsp_id for s in stats_list])
    _hline()
    print("Without Spread Factor")
    _num("NMD",           [s.nmd       for s in stats_list])
    _num("NMD10",         [s.nmd10     for s in stats_list])
    _num("NMD90",         [s.nmd90     for s in stats_list])
    _num("VMD",           [s.vmd       for s in stats_list])
    _num("VMD10",         [s.vmd10     for s in stats_list])
    _num("VMD90",         [s.vmd90     for s in stats_list])
    _num("CH",            [s.ch        for s in stats_list])
    _num("Mean diameter", [s.mean_diam for s in stats_list])
    _num("Diameter STD",  [s.std_diam  for s in stats_list])
    _hline()
    print("With Spread Factor")
    _num("NMD",           [s.nmd_sf       for s in stats_list])
    _num("NMD10",         [s.nmd10_sf     for s in stats_list])
    _num("NMD90",         [s.nmd90_sf     for s in stats_list])
    _num("VMD",           [s.vmd_sf       for s in stats_list])
    _num("VMD10",         [s.vmd10_sf     for s in stats_list])
    _num("VMD90",         [s.vmd90_sf     for s in stats_list])
    _num("CH",            [s.ch_sf        for s in stats_list])
    _num("Mean diameter", [s.mean_diam_sf for s in stats_list])
    _num("Diameter STD",  [s.std_diam_sf  for s in stats_list])
    _hline()
    _pct("Covered area",      [s.covered_area_pct for s in stats_list])
    _int("Droplets per card", [s.droplet_count    for s in stats_list])
    _num("Droplets per cm2",  [s.droplets_per_cm2 for s in stats_list])
    _hline()
    print()


# ── Per-card plot data ────────────────────────────────────────────────────────

@dataclass
class PlotData:
    stain_diam:    np.ndarray   # stain diameter per area class (µm)
    sf_diam:       np.ndarray   # spread-factor corrected diameter per area class (µm)
    area_hist:     np.ndarray   # droplet count per area class
    csn:           np.ndarray   # cumulative number % vs stain_diam
    bin_labels:    list          # e.g. ["<50", "50–100", …, ">600"]
    bin_counts:    np.ndarray   # droplet counts per bin (raw diameters)
    bin_counts_sf: np.ndarray   # droplet counts per bin (SF corrected diameters)


def compute_plot_data(
    droplets: list[Droplet],
    dpi: int,
    cfg: Config,
) -> PlotData | None:
    """Compute arrays needed to draw histogram and cumulative figures.

    Returns ``None`` when no droplets are detected.
    """
    if not droplets:
        return None

    pixel_um = 25_400.0 / dpi
    areas    = np.array([d.area_px for d in droplets], dtype=int)
    max_area = int(areas.max())
    area_hist = np.bincount(areas, minlength=max_area + 1)[1:]

    idx        = np.arange(1, len(area_hist) + 1, dtype=float)
    stain_diam = pixel_um * np.sqrt(4.0 * idx / np.pi)
    sf_diam    = apply_spread_factor(stain_diam)

    csn = 100.0 * np.cumsum(area_hist) / area_hist.sum()

    # Bin labels matching MATLAB's limiti_diametri convention
    edges      = cfg.hist_bin_edges
    bin_labels = [f"<{edges[0]:.0f}"]
    for i in range(len(edges) - 1):
        bin_labels.append(f"{edges[i]:.0f}–{edges[i+1]:.0f}")
    bin_labels.append(f">{edges[-1]:.0f}")

    def _bin(diams: np.ndarray) -> np.ndarray:
        counts    = np.zeros(len(bin_labels), dtype=float)
        counts[0] = area_hist[diams < edges[0]].sum()
        for i in range(len(edges) - 1):
            mask       = (diams >= edges[i]) & (diams < edges[i + 1])
            counts[i + 1] = area_hist[mask].sum()
        counts[-1] = area_hist[diams >= edges[-1]].sum()
        return counts

    return PlotData(
        stain_diam    = stain_diam,
        sf_diam       = sf_diam,
        area_hist     = area_hist,
        csn           = csn,
        bin_labels    = bin_labels,
        bin_counts    = _bin(stain_diam),
        bin_counts_sf = _bin(sf_diam),
    )


def save_scan_plots(
    plot_data_list: list[PlotData | None],
    stats_list: list[CardStats],
    out_path: Path,
    cfg: Config,
) -> None:
    """Save four figures for a whole scan — raw and SF-corrected kept separate.

    hist_raw_<stem>.png  : grouped bar chart, raw stain diameters, one bar per card
    hist_sf_<stem>.png   : grouped bar chart, SF-corrected drop diameters
    cumul_raw_<stem>.png : cumulative count % vs stain diameter, one line per card
    cumul_sf_<stem>.png  : cumulative count % vs SF-corrected drop diameter
    """
    # MATLAB colour palette (one colour per card)
    _COLOURS = [
        "#F44336", "#2196F3", "#4CAF50", "#00BCD4",
        "#9C27B0", "#FF9800", "#607D8B", "#E91E63",
        "#CDDC39", "#FF5722", "#3F51B5", "#009688",
        "#795548",
    ]

    valid = [(pd, st) for pd, st in zip(plot_data_list, stats_list) if pd is not None]
    if not valid:
        return

    stem    = out_path.stem   # e.g. "plots_card_example"
    out_dir = out_path.parent

    bin_labels = valid[0][0].bin_labels
    n_bins     = len(bin_labels)
    n_cards    = len(valid)
    bar_w      = 0.8 / n_cards

    def _save_hist(counts_key: str, title: str, out_file: Path) -> None:
        fig, ax = plt.subplots(figsize=(max(10, n_bins * 1.5), 8))
        fig.suptitle(title, fontsize=12)
        x = np.arange(n_bins)
        for i, (pd_i, st) in enumerate(valid):
            colour = _COLOURS[i % len(_COLOURS)]
            offset = (i - n_cards / 2 + 0.5) * bar_w
            counts = pd_i.bin_counts if counts_key == "raw" else pd_i.bin_counts_sf
            ax.bar(x + offset, counts, width=bar_w * 0.9,
                   color=colour, alpha=0.85, label=st.wsp_id,
                   edgecolor="black", linewidth=0.5)
        ax.set_xticks(x)
        ax.set_xticklabels(bin_labels, rotation=30, ha="right", fontsize=8)
        ax.set_xlabel("Diameter (µm)")
        ax.set_ylabel("Droplet count")
        ax.legend(fontsize=8, loc="upper right")
        plt.tight_layout()
        fig.savefig(out_file, dpi=cfg.dpi)
        plt.close(fig)
        print(f"  Histogram saved -> {out_file}")

    def _save_cumul(diam_key: str, title: str, out_file: Path) -> None:
        max_d = cfg.cumulative_max_diam_um
        fig, ax = plt.subplots(figsize=(10, 8))
        fig.suptitle(title, fontsize=12)
        for i, (pd_i, st) in enumerate(valid):
            colour = _COLOURS[i % len(_COLOURS)]
            diams  = pd_i.stain_diam if diam_key == "raw" else pd_i.sf_diam
            mask   = diams <= max_d
            ax.plot(
                np.concatenate([[0], diams[mask]]),
                np.concatenate([[0], pd_i.csn[mask]]),
                color=colour, linewidth=1.2, label=st.wsp_id,
            )
        ax.set_xlim(0, max_d)
        ax.set_ylim(0, 100)
        ax.set_xlabel("Diameter (µm)")
        ax.set_ylabel("Cumulative count (%)")
        ax.legend(fontsize=8, loc="lower right")
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        fig.savefig(out_file, dpi=cfg.dpi)
        plt.close(fig)
        print(f"  Cumulative saved -> {out_file}")

    _save_hist("raw", f"Droplets Count (stain diameter) — {stem}",
               out_dir / f"hist_raw_{stem}.png")
    _save_hist("sf",  f"Droplets Count (drop diameter, SF) — {stem}",
               out_dir / f"hist_sf_{stem}.png")
    _save_cumul("raw", f"Cumulative Sum (stain diameter) — {stem}",
                out_dir / f"cumul_raw_{stem}.png")
    _save_cumul("sf",  f"Cumulative Sum (drop diameter, SF) — {stem}",
                out_dir / f"cumul_sf_{stem}.png")


