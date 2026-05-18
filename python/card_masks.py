"""card_masks.py — Stage 1: WSP card segmentation.

Reads a scan image, detects each Water Sensitive Paper (WSP) card, and saves
the result as an overlay image showing each card border drawn in red.

The caller converts RGB → LAB once (in LAB mode) or skips conversion (in RGB
mode) and passes the working image here; the RGB image is only needed for
visualisation.

Two colour-space modes are supported (set via ``colour_space`` in config.yaml):

  ``lab`` (default)
      Perceptual LAB space.  Yellow paper detected by b* > yellow_b_min;
      blue stains by b* < blue_b_max AND chroma > min_chroma.  The chroma
      guard rejects dark achromatic inter-card shadows.

  ``rgb``
      Direct RGB channels, matching the original MATLAB pipeline.
      Yellow paper: R > B + yellow_rg_over_b  AND  G > B + yellow_rg_over_b.
      Blue stains : B > R + blue_b_over_rg    AND  B > G + blue_b_over_rg.

Long thin edge streaks (raised-card shadows, humidity) are handled in Stage 2
by eroding the card mask before droplet detection.

Configuration
-------------
All parameters are read from ``config.yaml`` (one file, one load for the whole
pipeline).  Any value can be overridden via CLI arguments.

Usage
-----
    python card_masks.py <image_file> [--output <dir>] [--config config.yaml]
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import yaml
import tifffile
from PIL import Image
from skimage.color import rgb2lab
from skimage.morphology import binary_closing, binary_erosion, disk, remove_small_objects
from skimage.measure import label, regionprops, find_contours
from scipy.ndimage import binary_fill_holes


_DEFAULT_CONFIG = Path(__file__).parent / "config.yaml"


# ── Configuration ─────────────────────────────────────────────────────────────

@dataclass
class Config:
    # scanner
    dpi: int
    # colour space: 'lab' or 'rgb'
    colour_space: str
    # card detection — LAB thresholds
    yellow_b_min:     float
    blue_b_max:       float
    min_chroma:       float
    # card detection — RGB thresholds (MATLAB-compatible)
    yellow_rg_over_b: float
    blue_b_over_rg:   float
    # morphological cleanup
    min_card_area_px: int
    closing_radius:   int
    closing_se:       str   # 'disk' or 'matlab_disk'
    # droplet detection — LAB thresholds
    droplet_b_max:            float
    droplet_l_max:            float
    # droplet detection — RGB threshold (MATLAB-compatible)
    droplet_r_max:            float
    # droplet detection — shared settings
    min_droplet_area_px:      int
    subseg_area_threshold_px: int
    border_erosion_px:        int
    overlay_alpha:            float
    hist_bin_edges:           list
    cumulative_max_diam_um:   float
    # sub-segmentation
    subseg_method:                str
    adaptive_initial_threshold:   float
    adaptive_fraction:            float


def load_config(path: Path | None = None) -> Config:
    """Load the full pipeline config from a YAML file (called once)."""
    yaml_path = path or _DEFAULT_CONFIG
    with open(yaml_path) as f:
        raw = yaml.safe_load(f)

    s  = raw["scanner"]
    cs = str(raw.get("colour_space", "lab")).lower()
    if cs not in ("lab", "rgb"):
        raise ValueError(f"colour_space must be 'lab' or 'rgb', got: {cs!r}")
    c  = raw["card_detection"]
    m  = raw["morphology"]
    d  = raw["droplet_detection"]

    subseg_method = str(d.get("subseg_method", "watershed")).lower()
    if subseg_method not in ("watershed", "adaptive_threshold"):
        raise ValueError(
            f"subseg_method must be 'watershed' or 'adaptive_threshold', got: {subseg_method!r}"
        )
    adaptive_fraction = float(d.get("adaptive_fraction", 2))
    if adaptive_fraction <= 1:
        raise ValueError(
            f"adaptive_fraction must be > 1 (got {adaptive_fraction}); "
            "values ≤ 1 would push the threshold outside the intended range"
        )

    closing_se = str(m.get("closing_se", "disk")).lower()
    if closing_se not in ("disk", "matlab_disk"):
        raise ValueError(
            f"closing_se must be 'disk' or 'matlab_disk', got: {closing_se!r}"
        )

    return Config(
        dpi                      = int(s["dpi"]),
        colour_space             = cs,
        yellow_b_min             = float(c["yellow_b_min"]),
        blue_b_max               = float(c["blue_b_max"]),
        min_chroma               = float(c["min_chroma"]),
        yellow_rg_over_b         = float(c["yellow_rg_over_b"]),
        blue_b_over_rg           = float(c["blue_b_over_rg"]),
        min_card_area_px         = int(m["min_card_area_px"]),
        closing_radius           = int(m["closing_radius"]),
        closing_se               = closing_se,
        droplet_b_max            = float(d["droplet_b_max"]),
        droplet_l_max            = float(d["droplet_l_max"]),
        droplet_r_max            = float(d["droplet_r_max"]),
        min_droplet_area_px      = int(d["min_droplet_area_px"]),
        subseg_area_threshold_px = int(d["subseg_area_threshold_px"]),
        border_erosion_px        = int(d["border_erosion_px"]),
        overlay_alpha            = float(d["overlay_alpha"]),
        hist_bin_edges           = [float(v) for v in d["hist_bin_edges"]],
        cumulative_max_diam_um   = float(d["cumulative_max_diam_um"]),
        subseg_method              = subseg_method,
        adaptive_initial_threshold = float(d.get("adaptive_initial_threshold", 190)),
        adaptive_fraction          = adaptive_fraction,
    )


# ── Data model ────────────────────────────────────────────────────────────────

@dataclass
class CardMask:
    label:    str           # e.g. "WSP-1"
    mask:     np.ndarray    # bool H×W — eroded card area (used for all analysis)
    bbox:     tuple         # (row_min, col_min, row_max, col_max)
    centroid: tuple         # (row, col)


# ── Public API ────────────────────────────────────────────────────────────────

def create_card_masks(
    img: np.ndarray,
    card_names: list[str] | None = None,
    cfg: Config | None = None,
) -> list[CardMask]:
    """Segment WSP cards from the working image.

    Parameters
    ----------
    img:
        H×W×3 working image.  In LAB mode (``cfg.colour_space == 'lab'``) this
        must be a float32 LAB image; in RGB mode it must be a uint8 RGB image.
        The caller performs any conversion once before calling this function.
    card_names:
        Optional custom labels; generic labels used for any missing entry.
    cfg:
        Pipeline configuration.

    Returns
    -------
    list[CardMask]
        Ordered left-column top->bottom, then right-column top->bottom.
    """
    raw_mask   = _build_card_mask(img, cfg)
    clean_mask = _cleanup(raw_mask, cfg)
    regions    = _get_regions(clean_mask)
    return _order_cards(regions, clean_mask, card_names, cfg)


# ── Internal helpers ──────────────────────────────────────────────────────────

def _build_card_mask(img: np.ndarray, cfg: Config) -> np.ndarray:
    """Binary mask covering WSP card regions; dispatches to LAB or RGB mode."""
    if cfg.colour_space == "rgb":
        return _build_card_mask_rgb(img, cfg)
    return _build_card_mask_lab(img, cfg)


def _build_card_mask_lab(lab: np.ndarray, cfg: Config) -> np.ndarray:
    """Binary mask from LAB channels: yellow paper OR chromatic blue stains."""
    b      = lab[:, :, 2]
    a      = lab[:, :, 1]
    chroma = np.sqrt(a**2 + b**2)

    yellow = b > cfg.yellow_b_min
    # Chroma guard: achromatic inter-card shadows have chroma ~0-5.
    blue   = (b < cfg.blue_b_max) & (chroma > cfg.min_chroma)

    return yellow | blue


def _build_card_mask_rgb(rgb: np.ndarray, cfg: Config) -> np.ndarray:
    """Binary mask from RGB channels: yellow paper OR blue stains.

    Matches MATLAB ``creaMascheraCartina_n``:
      yellow = R > B + yellow_rg_over_b  AND  G > B + yellow_rg_over_b
      blue   = B > R + blue_b_over_rg    AND  B > G + blue_b_over_rg
    """
    r = rgb[:, :, 0].astype(np.int16)
    g = rgb[:, :, 1].astype(np.int16)
    b = rgb[:, :, 2].astype(np.int16)

    yellow = (r > b + cfg.yellow_rg_over_b) & (g > b + cfg.yellow_rg_over_b)
    blue   = (b > r + cfg.blue_b_over_rg)   & (b > g + cfg.blue_b_over_rg)

    return yellow | blue


def _matlab_strel_disk(r: int, n: int = 8) -> np.ndarray:
    """Return a structuring element matching MATLAB's ``strel('disk', r)``.

    MATLAB's default (n=8) decomposes the disk into eight line SEs and takes
    their union, producing a 16-sided convex polygon that circumscribes the
    true circle.  At the 22.5-degree bisector angles between adjacent lines the
    polygon radius is ``r / cos(π/16) ≈ r * 1.0196`` — roughly one extra pixel
    for r=50.  The true-disk SE (skimage ``disk(r)``) is therefore slightly
    *smaller* than this polygon, causing a 2–4 droplet difference when the extra
    border strip captures additional stain pixels.
    """
    size = 2 * r + 5
    c = r + 2
    yy, xx = np.mgrid[0:size, 0:size]
    yy = yy - c
    xx = xx - c
    se = np.ones((size, size), dtype=bool)
    for k in range(n):
        angle = k * np.pi / n
        proj = xx * np.cos(angle) + yy * np.sin(angle)
        se &= (np.abs(proj) <= r)
    return se


def _cleanup(mask: np.ndarray, cfg: Config) -> np.ndarray:
    """Remove noise, close gaps, and fill internal holes.

    Order matches MATLAB (creaMascheraCartina_n):
      1. imfill(mask, 'holes')      → binary_fill_holes
      2. bwareaopen(mask, min_size) → remove_small_objects
      3. imclose(mask, strel(r))    → binary_closing

    The structuring element shape is controlled by ``cfg.closing_se``:
      * ``'disk'``        — true Euclidean disk (skimage ``disk(r)``); default
      * ``'matlab_disk'`` — 16-sided polygon matching MATLAB ``strel('disk', r)``
    """
    mask = binary_fill_holes(mask)
    mask = remove_small_objects(mask, min_size=cfg.min_card_area_px)
    se = (
        _matlab_strel_disk(cfg.closing_radius)
        if cfg.closing_se == "matlab_disk"
        else disk(cfg.closing_radius)
    )
    mask = binary_closing(mask, se)
    return mask


def _erode_mask(mask: np.ndarray, cfg: Config) -> np.ndarray:
    """Erode each connected card region to exclude border artefacts."""
    return binary_erosion(mask, disk(cfg.border_erosion_px))


def _get_regions(mask: np.ndarray):
    labeled = label(mask, connectivity=2)
    return regionprops(labeled)


def _count_axis_clusters(values: np.ndarray, dimension: int,
                         gap_fraction: float = 0.10) -> int:
    """Return the number of clusters along one axis using a relative gap threshold.

    A gap between two consecutive sorted values is treated as a cluster boundary
    when it exceeds ``gap_fraction * dimension``.
    """
    sorted_vals = np.sort(values.astype(float))
    if len(sorted_vals) <= 1:
        return 1
    gaps = np.diff(sorted_vals)
    threshold = gap_fraction * dimension
    return 1 + int(np.sum(gaps > threshold))


def _assign_clusters(values: np.ndarray, dimension: int,
                     gap_fraction: float = 0.10) -> np.ndarray:
    """Return a cluster index (0-based) for each element of *values*.

    Clusters are ordered by ascending mean value (top→bottom or left→right).
    """
    order   = np.argsort(values.astype(float))
    labels  = np.empty(len(values), dtype=int)
    threshold = gap_fraction * dimension
    cluster = 0
    labels[order[0]] = 0
    for k in range(1, len(order)):
        if values[order[k]] - values[order[k - 1]] > threshold:
            cluster += 1
        labels[order[k]] = cluster
    return labels


def _detect_layout(centroids: list[tuple], img_h: int, img_w: int) -> str:
    """Detect how cards are arranged and return the layout name.

    Returns
    -------
    'rows'    — 1 or 2 horizontal rows detected (order: left→right, top→bottom)
    'columns' — 1 or 2 vertical columns detected (order: top→bottom, left→right)
    'other'   — neither: fall back to centroid top→bottom sort

    When both axes independently show 1–2 clusters (e.g. a 2×2 grid), rows
    take priority.
    """
    ys = np.array([c[0] for c in centroids])
    xs = np.array([c[1] for c in centroids])

    n_rows = _count_axis_clusters(ys, img_h)
    n_cols = _count_axis_clusters(xs, img_w)

    row_ok = n_rows in (1, 2)
    col_ok = n_cols in (1, 2)

    if row_ok:           # rows take priority (also covers both-detected case)
        return "rows"
    elif col_ok:
        return "columns"
    else:
        return "other"


def _order_cards(regions, mask: np.ndarray, card_names, cfg: Config) -> list[CardMask]:
    """Detect card layout and order cards accordingly, then assign labels.

    Layout detection
    ----------------
    rows    (1 or 2 y-clusters): order left→right within each row, top→bottom
    columns (1 or 2 x-clusters): order top→bottom within each column, left→right
    other                       : order by y-centroid top→bottom
    """
    if not regions:
        return []

    H, W    = mask.shape
    labeled = label(mask, connectivity=2)

    centroids = [r.centroid for r in regions]
    layout    = _detect_layout(centroids, H, W)

    _layout_labels = {
        "rows":    f"row layout ({_count_axis_clusters(np.array([c[0] for c in centroids]), H)} row(s),"
                   f" ordered left-to-right then top-to-bottom)",
        "columns": f"column layout ({_count_axis_clusters(np.array([c[1] for c in centroids]), W)} column(s),"
                   f" ordered top-to-bottom then left-to-right)",
        "other":   "no row/column structure detected, ordered top-to-bottom by centroid",
    }
    print(f"  Layout detected: {_layout_labels[layout]}")

    if layout == "rows":
        ys      = np.array([c[0] for c in centroids])
        row_idx = _assign_clusters(ys, H)
        ordered = [r for _, r in sorted(enumerate(regions),
                                        key=lambda ir: (row_idx[ir[0]],
                                                        ir[1].centroid[1]))]
    elif layout == "columns":
        xs      = np.array([c[1] for c in centroids])
        col_idx = _assign_clusters(xs, W)
        ordered = [r for _, r in sorted(enumerate(regions),
                                        key=lambda ir: (col_idx[ir[0]],
                                                        ir[1].centroid[0]))]
    else:
        ordered = sorted(regions, key=lambda r: r.centroid[0])

    result = []
    for i, r in enumerate(ordered):
        lbl       = (card_names[i] if card_names and i < len(card_names)
                     else f"WSP-{i + 1}")
        full_mask = (labeled == r.label)
        eroded    = _erode_mask(full_mask, cfg)
        result.append(CardMask(
            label    = lbl,
            mask     = eroded,
            bbox     = r.bbox,
            centroid = r.centroid,
        ))
    return result


# ── Image I/O ─────────────────────────────────────────────────────────────────

def load_image(path: Path) -> np.ndarray:
    """Load a scan image file as a uint8 RGB array.

    16-bit images are scaled by the dtype maximum (65535) so that intensity
    values are consistent across images and RGB thresholds remain stable.
    """
    suffix = path.suffix.lower()
    if suffix in (".tif", ".tiff"):
        img = tifffile.imread(str(path))
    else:
        img = np.array(Image.open(path))

    if img.dtype != np.uint8:
        max_val = np.iinfo(img.dtype).max if np.issubdtype(img.dtype, np.integer) else img.max()
        img = (img.astype(np.float32) / max_val * 255).astype(np.uint8)
    if img.ndim == 2:
        img = np.stack([img] * 3, axis=-1)
    if img.shape[2] == 4:
        img = img[:, :, :3]
    return img


# ── Visualisation ─────────────────────────────────────────────────────────────

def _draw_contours(ax, masks: list[np.ndarray], color="red", linewidth=1):
    """Draw the boundary contour of each boolean mask onto *ax*."""
    for mask in masks:
        for contour in find_contours(mask, level=0.5):
            ax.plot(contour[:, 1], contour[:, 0],
                    color=color, linewidth=linewidth)


def _draw_mask_overlay(ax, masks: list[np.ndarray], alpha: float = 0.5) -> None:
    """Overlay a semitransparent red fill over the *background* (everything
    outside the union of *masks*), leaving the card regions clearly visible."""
    if not masks:
        return
    H, W = masks[0].shape
    combined = np.zeros((H, W), dtype=bool)
    for m in masks:
        combined |= m
    rgba = np.zeros((H, W, 4), dtype=np.float32)
    rgba[~combined] = [1.0, 0.0, 1.0, alpha]   # background → red
    ax.imshow(rgba)


def save_overlay(rgb: np.ndarray, cards: list[CardMask], out_path: Path, cfg: Config) -> None:
    """Save the scan with each card border drawn in red."""
    fig, ax = plt.subplots(figsize=(12, 16))
    ax.imshow(rgb)
    ax.axis("off")
    ax.set_title(f"Card masks — {out_path.stem}  ({len(cards)} cards found)")

    _draw_contours(ax, [c.mask for c in cards])

    for card in cards:
        cy, cx = card.centroid
        ax.text(cx, cy, card.label,
                color="white", fontsize=11, ha="center", va="center",
                fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.3", facecolor="red",
                          alpha=0.7, edgecolor="none"))

    plt.tight_layout()
    fig.savefig(out_path, dpi=cfg.dpi)
    plt.close(fig)
    print(f"  Overlay saved -> {out_path}")


def save_mask_stages(
    rgb: np.ndarray,
    raw_mask: np.ndarray,
    clean_mask: np.ndarray,
    cards: list[CardMask],
    out_path: Path,
    cfg: Config,
) -> None:
    """Save a 1x3 figure: raw LAB threshold -> cleanup -> ordered cards."""
    labeled_clean = label(clean_mask, connectivity=2)
    region_masks  = [labeled_clean == r.label for r in regionprops(labeled_clean)]

    fig, axes = plt.subplots(1, 3, figsize=(24, 10))
    fig.suptitle(f"Mask pipeline — {out_path.stem}", fontsize=14)

    axes[0].imshow(raw_mask, cmap="gray")
    threshold_label = (
        "LAB threshold\n(yellow | chromatic blue)"
        if cfg.colour_space == "lab"
        else "RGB threshold\n(yellow | blue)"
    )
    axes[0].set_title(f"1. {threshold_label}")
    axes[0].axis("off")

    axes[1].imshow(rgb)
    _draw_mask_overlay(axes[1], region_masks)
    axes[1].set_title(f"2. After cleanup\n({len(region_masks)} regions)")
    axes[1].axis("off")

    axes[2].imshow(rgb)
    _draw_mask_overlay(axes[2], [c.mask for c in cards])
    axes[2].set_title(f"3. Ordered cards (eroded)\n({len(cards)} cards)")
    axes[2].axis("off")
    for card in cards:
        cy, cx = card.centroid
        axes[2].text(cx, cy, card.label,
                     color="white", fontsize=9, ha="center", va="center",
                     fontweight="bold",
                     bbox=dict(boxstyle="round,pad=0.2", facecolor="red",
                               alpha=0.7, edgecolor="none"))

    plt.tight_layout()
    fig.savefig(out_path, dpi=cfg.dpi)
    plt.close(fig)
    print(f"  Stage figure saved -> {out_path}")


def save_lab_bands(lab: np.ndarray, out_path: Path, cfg: Config) -> None:
    """Save a 1x3 figure showing the L*, a*, b* channels of the LAB image."""
    titles = ["L* (lightness)", "a* (green<->red)", "b* (blue<->yellow)"]
    cmaps  = ["gray", "gray", "gray"]

    fig, axes = plt.subplots(1, 3, figsize=(24, 10))
    fig.suptitle(f"LAB channels — {out_path.stem}\n", fontsize=14)

    for i, (title, cmap) in enumerate(zip(titles, cmaps)):
        im = axes[i].imshow(lab[:, :, i], cmap=cmap)
        axes[i].set_title(title)
        axes[i].axis("off")
        plt.colorbar(im, ax=axes[i], fraction=0.046, pad=0.04)

    plt.tight_layout()
    fig.savefig(out_path, dpi=cfg.dpi)
    plt.close(fig)
    print(f"  LAB bands saved -> {out_path}")


def save_rgb_channels(rgb: np.ndarray, out_path: Path, cfg: Config) -> None:
    """Save a 1x3 figure showing the R, G, B channels of the RGB image."""
    titles = ["R (red)", "G (green)", "B (blue)"]
    cmaps  = ["Reds_r", "Greens_r", "Blues_r"]

    fig, axes = plt.subplots(1, 3, figsize=(24, 10))
    fig.suptitle(f"RGB channels — {out_path.stem}\n", fontsize=14)

    for i, (title, cmap) in enumerate(zip(titles, cmaps)):
        im = axes[i].imshow(rgb[:, :, i], cmap=cmap, vmin=0, vmax=255)
        axes[i].set_title(title)
        axes[i].axis("off")
        plt.colorbar(im, ax=axes[i], fraction=0.046, pad=0.04)

    plt.tight_layout()
    fig.savefig(out_path, dpi=cfg.dpi)
    plt.close(fig)
    print(f"  RGB channels saved -> {out_path}")


