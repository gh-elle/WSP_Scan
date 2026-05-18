"""detect_droplets.py — Stage 2: droplet detection on individual WSP cards.

For each card mask produced by Stage 1, this module:
  1. Crops the full-scan working image to the card bounding box.
  2. Uses the eroded card mask (already applied in Stage 1).
  3. Binarises using the active colour-space mode:
       LAB mode: b* < droplet_b_max AND L* < droplet_l_max
       RGB mode: R < droplet_r_max  (matches MATLAB binarizza_cartina3)
  4. Finds connected components; discards noise below min_droplet_area_px.
  5. Sub-segments large blobs (probable merged drops) with watershed.

The caller converts RGB -> LAB once (in LAB mode) or skips conversion (in RGB
mode) and passes the working image here; the RGB image is passed only to the
visualisation functions.

Configuration
-------------
Parameters are read from the ``droplet_detection`` section of ``config.yaml``
via the shared ``load_config`` / ``Config`` from card_masks.py.

Usage
-----
    python detect_droplets.py <image_file> [--output <dir>] [--config config.yaml]
"""

from __future__ import annotations

from collections import namedtuple
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from skimage.color import rgb2lab
from skimage.measure import label, regionprops
from skimage.segmentation import watershed
from skimage.feature import peak_local_max
from skimage.morphology import disk
from scipy.ndimage import distance_transform_edt

from card_masks import (
    CardMask, Config, load_config,
    _build_card_mask, _cleanup, _get_regions, _order_cards,
    load_image,
)


# ── Data model ────────────────────────────────────────────────────────────────

# Lightweight stand-in for a regionprops region; used when sub-segmentation
# returns components in a different coordinate frame (adaptive threshold mode).
_SubRegion = namedtuple("_SubRegion", ["area", "centroid"])


@dataclass
class Droplet:
    area_px:  int    # stain area in pixels
    centroid: tuple  # (row, col) in card-crop coordinates


# ── Public API ────────────────────────────────────────────────────────────────

def detect_droplets(
    img: np.ndarray,
    card: CardMask,
    cfg: Config,
) -> list[Droplet]:
    """Detect water-drop stains on a single WSP card.

    Parameters
    ----------
    img:
        Full scan as H×W×3 working image.  In LAB mode this is a float32 LAB
        array; in RGB mode it is a uint8 RGB array.
    card:
        Card mask from Stage 1 (defines the region of interest).
    cfg:
        Pipeline configuration (shared Config from card_masks).

    Returns
    -------
    tuple[list[Droplet], int]
        Detected droplets after sub-segmentation, and the raw stain pixel count
        (used for MATLAB-compatible coverage: ``R < 120`` pixels before subseg).
    """
    r0, c0, r1, c1 = card.bbox
    img_crop  = img[r0:r1, c0:c1]
    mask_crop = card.mask[r0:r1, c0:c1]   # already eroded in Stage 1

    stain_map  = _build_stain_map(img_crop, mask_crop, cfg)
    stain_px   = int(stain_map.sum())   # raw stain pixel count (for MATLAB-compatible coverage)

    components = _label_components(stain_map, cfg.min_droplet_area_px)
    components = _subsegment(stain_map, img_crop, components, cfg)

    droplets = [Droplet(area_px=int(r.area), centroid=r.centroid) for r in components]
    return droplets, stain_px


# ── Internal helpers ──────────────────────────────────────────────────────────

def _build_stain_map(
    img_crop: np.ndarray,
    eroded_mask: np.ndarray,
    cfg: Config,
) -> np.ndarray:
    """Boolean stain map inside the eroded card mask; dispatches to LAB or RGB mode."""
    if cfg.colour_space == "rgb":
        return _build_stain_map_rgb(img_crop, eroded_mask, cfg)
    return _build_stain_map_lab(img_crop, eroded_mask, cfg)


def _build_stain_map_lab(
    lab_crop: np.ndarray,
    eroded_mask: np.ndarray,
    cfg: Config,
) -> np.ndarray:
    """Stain map using LAB channels: blue-shifted, non-bright pixels inside mask.

    Conditions:
    - b* < droplet_b_max   : pixel is blue-shifted (primary colour selector)
    - L* < droplet_l_max   : pixel is not too bright (excludes white scratches)
    """
    return (
        (lab_crop[:, :, 2] < cfg.droplet_b_max) &
        (lab_crop[:, :, 0] < cfg.droplet_l_max) &
        eroded_mask
    )


def _build_stain_map_rgb(
    rgb_crop: np.ndarray,
    eroded_mask: np.ndarray,
    cfg: Config,
) -> np.ndarray:
    """Stain map using RGB channels: R < droplet_r_max inside mask.

    Matches MATLAB ``binarizza_cartina3``: ``bw = (R < soglia) & maschera``.
    """
    return (rgb_crop[:, :, 0] < cfg.droplet_r_max) & eroded_mask


def _label_components(bw: np.ndarray, min_area: int):
    """Label connected components; discard noise below min_area."""
    lbl = label(bw, connectivity=2)
    return [r for r in regionprops(lbl) if r.area >= min_area]


def _subsegment(stain_map: np.ndarray, img_crop: np.ndarray, regions, cfg: Config):
    """Split large blobs (probable merged drops) using the configured method.

    Small components are returned unchanged.  For components above
    subseg_area_threshold_px the selected method is applied:

    'watershed'           — distance-transform watershed (both LAB and RGB)
    'adaptive_threshold'  — MATLAB-compatible adaptive channel threshold
                            (RGB mode only; falls back to watershed in LAB mode)
    """
    result  = []
    lbl_map = label(stain_map, connectivity=2)

    for r in regions:
        use_adaptive = (
            cfg.subseg_method == "adaptive_threshold"
            and cfg.colour_space == "rgb"
        )
        if use_adaptive:
            # MATLAB isola_elementi3 processes ALL components regardless of size
            result.extend(_subseg_adaptive(img_crop, lbl_map, r, cfg))
        elif r.area < cfg.subseg_area_threshold_px:
            result.append(r)
        else:
            result.extend(_subseg_watershed(lbl_map, r, cfg))

    return result


def _subseg_watershed(lbl_map: np.ndarray, r, cfg: Config) -> list:
    """Watershed sub-segmentation on one component."""
    component_mask = lbl_map == r.label
    dist           = distance_transform_edt(component_mask)

    min_dist = max(3, int(np.sqrt(r.area / np.pi) * 0.3))
    peaks    = peak_local_max(dist, min_distance=min_dist, labels=component_mask)

    if len(peaks) <= 1:
        return [r]

    markers = np.zeros_like(dist, dtype=int)
    for i, (pr, pc) in enumerate(peaks, start=1):
        markers[pr, pc] = i

    ws = watershed(-dist, markers, mask=component_mask)
    result = [r2 for r2 in regionprops(ws) if r2.area >= cfg.min_droplet_area_px]
    return result if result else [r]


def _subseg_adaptive(
    img_crop: np.ndarray,
    lbl_map: np.ndarray,
    r,
    cfg: Config,
) -> list:
    """MATLAB-compatible adaptive R-channel threshold sub-segmentation.

    Matches ``isola_elementi3`` logic:
      threshold = soglia_iniziale - (soglia_iniziale - min(R_in_bbox)) / frazione

    The R channel of the component's bounding box is thresholded (no component
    masking, matching MATLAB).  Resulting sub-region centroids are offset back
    to card-crop coordinates so that ``Droplet.centroid`` remains consistent.

    Only called in RGB mode.  In LAB mode the caller falls back to watershed.
    """
    rb0, cb0, rb1, cb1 = r.bbox
    r_bbox = img_crop[rb0:rb1, cb0:cb1, 0].astype(float)   # R channel crop

    min_val   = r_bbox.min()
    threshold = cfg.adaptive_initial_threshold - (
        cfg.adaptive_initial_threshold - min_val
    ) / cfg.adaptive_fraction

    sub_bw  = r_bbox < threshold
    sub_lbl = label(sub_bw, connectivity=2)
    sub_regions = [
        r2 for r2 in regionprops(sub_lbl)
        if r2.area >= cfg.min_droplet_area_px
    ]

    if not sub_regions:
        return [r]

    # Offset centroids from bbox-local → card-crop coordinates
    return [
        _SubRegion(
            area     = r2.area,
            centroid = (r2.centroid[0] + rb0, r2.centroid[1] + cb0),
        )
        for r2 in sub_regions
    ]


# ── Visualisation ─────────────────────────────────────────────────────────────

def save_card_detections(
    rgb: np.ndarray,
    img: np.ndarray,
    cards: list[CardMask],
    droplets_per_card: list[list[Droplet]],
    cfg: Config,
    out_path: Path,
) -> None:
    """Save one image per card with detected stains as a semi-transparent red overlay.

    The red overlay is built from the binary stain map that matches exactly what
    detect_droplets uses (LAB or RGB mode, as configured).  The droplet count
    after sub-segmentation is shown in the title.

    Files are named ``<out_path.stem>_<card_label>.png``.
    """
    out_dir = out_path.parent
    stem    = out_path.stem

    for card, droplets in zip(cards, droplets_per_card):
        r0, c0, r1, c1 = card.bbox
        rgb_crop  = rgb[r0:r1, c0:c1]
        img_crop  = img[r0:r1, c0:c1]
        mask_crop = card.mask[r0:r1, c0:c1]   # already eroded in Stage 1
        stain_map = _build_stain_map(img_crop, mask_crop, cfg)

        # Build an RGBA overlay: red where stain, transparent elsewhere
        overlay = np.zeros((*stain_map.shape, 4), dtype=np.float32)
        overlay[stain_map, 0] = 1.0   # R
        overlay[stain_map, 3] = cfg.overlay_alpha

        fig, ax = plt.subplots(figsize=(8, 14))
        ax.imshow(rgb_crop)
        ax.imshow(overlay)
        ax.set_title(f"{card.label} — {len(droplets)} droplets")
        ax.axis("off")

        plt.tight_layout()
        safe_label = card.label.replace(" ", "_")
        print(f"Saving card detection: {stem}_{safe_label}.png ...")
        fig.savefig(out_dir / f"{stem}_{safe_label}.png", dpi=cfg.dpi)
        plt.close(fig)

    print(f"  Card detection images saved -> {out_dir}/{stem}_*.png")


