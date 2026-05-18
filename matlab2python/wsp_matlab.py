"""wsp_matlab.py — Direct Python translation of the MATLAB WSP pipeline.

Mirrors the MATLAB pipeline in matlab/:
    scan_cards.m → analyze_card_image.m → compute_droplet_coverage.m

The algorithm, thresholds, and constants are identical to those in the MATLAB code.
Only external parameter: DPI (default 600).  No config file.

Usage
-----
    python wsp_matlab.py <image.tif>                               # stdout only
    python wsp_matlab.py <image.tif> --dpi 1200
    python wsp_matlab.py <image.tif> --results-folder ./results    # saves CSV + PNGs
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import tifffile
from scipy.ndimage import binary_fill_holes
from skimage.measure import label, regionprops
from skimage.morphology import binary_closing, disk, remove_small_objects

# ── Constants ─────────────────────────────────────────────────────────────────

CARD_AREA_CM2        = 2.6 * 7.6   # droplet_statistics.m: card_area = 2.6 * 7.6
STAIN_THRESHOLD      = 120         # binarize_card.m: threshold = 120
ADAPTIVE_INITIAL     = 190         # isolate_elements.m: initial_threshold = 190
ADAPTIVE_FRACTION    = 2           # isolate_elements.m: fraction = 2
MIN_CARD_AREA        = 1000        # create_card_mask.m: bwareaopen(card_mask, 1000)
DISK_RADIUS          = 50          # create_card_mask.m: strel('disk', 50)
DIAMETER_LIMITS      = [50, 100, 150, 200, 300, 400, 500, 600]  # analyze_card_image.m
COLORS               = [                                          # analyze_card_image.m
    (1, 0, 0), (0, 0, 1), (0, 1, 0), (0, 1, 1), (1, 0, 1),
    (1, 0.5, 0), (0.5, 0.5, 1), (1, 0.5, 0.5), (1, 1, 0),
    (0.5, 1, 0.5), (0.3, 1, 0.3), (0.5, 0.5, 0.5),
]
CSV_HEADER           = [                                          # scan_cards.m columns
    'filename', 'card',
    'NMD', 'NMD10', 'NMD90', 'VMD', 'VMD10', 'VMD90', 'CH', 'mean_diameter', 'std_diameter',
    'covered_area_pct', 'droplet_count', 'droplets_per_cm2',
    'NMD_SF', 'NMD10_SF', 'NMD90_SF', 'VMD_SF', 'VMD10_SF', 'VMD90_SF', 'CH_SF',
    'mean_diameter_SF', 'std_diameter_SF',
]

# Path to the pre-computed exact MATLAB strel('disk', 50) kernel (optional).
# Generated with: a=zeros(200,200); a(100,100)=true; se=strel('disk',50);
#                 imwrite(uint8(imdilate(a,se))*255, 'matlab_strel_disk50.png')
# Then extracted with numpy and saved as .npy:
#   kernel = np.array(Image.open('matlab_strel_disk50.png'))[:,:,0] > 0
#   rows/cols = bounding box of True pixels → np.save('matlab_strel_disk50.npy', kernel)
# When present, this file is used instead of the Python approximation in
# build_card_mask(), giving results numerically identical to MATLAB.
_MATLAB_SE_NPY = Path(__file__).parent / "matlab_strel_disk50.npy"


# ── Stage 1: card mask (create_card_mask.m) ───────────────────────────────────

def _matlab_strel_disk(r: int, n: int = 4) -> np.ndarray:
    """Octagonal structuring element approximating MATLAB ``strel('disk', r, n)``.

    Used as fallback when the exact pre-computed MATLAB SE file is not available.
    MATLAB's strel('disk', r) uses N=4 periodic-line decomposition by default,
    producing an octagonal polygon via intersection of N stripe conditions at
    angles k*π/N for k=0…N-1.
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


def _load_card_se() -> np.ndarray:
    """Load the structuring element for card mask closing (create_card_mask.m).

    MATLAB's strel('disk', 50) with default N=4 decomposition produces a
    99×99 SE (effective radius 49) with rounded corners that extend to ~53 px.
    This differs from the simple octagon intersection formula used by
    ``_matlab_strel_disk``, which gives a 101×101 SE.

    When ``matlab_strel_disk50.npy`` is present alongside this script (generated
    by exporting MATLAB's actual strel to a PNG and converting with numpy), that
    exact kernel is used.  This gives per-card drop counts identical to MATLAB.
    When the file is absent the octagon approximation is used as fallback.
    """
    if _MATLAB_SE_NPY.exists():
        return np.load(str(_MATLAB_SE_NPY))
    return _matlab_strel_disk(DISK_RADIUS)


def build_card_mask(R: np.ndarray, G: np.ndarray, B: np.ndarray) -> np.ndarray:
    """Replicate create_card_mask.m mask logic.

    MATLAB (create_card_mask.m):
        card_mask = (R > B+50 & G > B+50) | (B > R+15 & B > G+15)
        card_mask = imfill(card_mask, 'holes')
        card_mask = bwareaopen(card_mask, 1000)
        card_mask = imclose(card_mask, strel('disk', 50))
    """
    R = R.astype(np.int16)
    G = G.astype(np.int16)
    B = B.astype(np.int16)

    mask = ((R > B + 50) & (G > B + 50)) | ((B > R + 15) & (B > G + 15))

    # imfill(mask, 'holes') — scipy default is fully connected (all directions)
    mask = binary_fill_holes(mask)

    # bwareaopen(mask, 1000) — MATLAB 2-D default is 8-connectivity
    mask = remove_small_objects(mask, min_size=MIN_CARD_AREA, connectivity=2)

    # imclose(mask, strel('disk', 50)) — uses exact MATLAB SE when available
    mask = binary_closing(mask, _load_card_se())

    return mask


def order_cards(regions, mask_shape: tuple) -> list:
    """Replicate card ordering from create_card_mask.m.

    MATLAB (create_card_mask.m):
        left_range = fix(size(card_mask, 2) / 3)
        left_idx   = find(x_pos <  left_range)
        right_idx  = find(x_pos >= left_range)
        [~, sort_left]  = sort(y_pos(left_idx))
        [~, sort_right] = sort(y_pos(right_idx))
        idx = [left_idx(sort_left), right_idx(sort_right)]

    x_pos = BoundingBox column coordinate (0-indexed in skimage)
    y_pos = BoundingBox row coordinate    (0-indexed in skimage)
    """
    # skimage bbox = (row_min, col_min, row_max, col_max) — 0-indexed
    x_pos = np.array([r.bbox[1] for r in regions], dtype=float)  # col
    y_pos = np.array([r.bbox[0] for r in regions], dtype=float)  # row

    left_range = int(mask_shape[1] / 3)   # fix() truncates toward zero

    # MATLAB's bwconncomp is column-major so left-column cards always come first
    # in the component list — `y_pos(1:n_left)` ARE the left cards.
    # skimage label() is row-major, so we must partition explicitly by x_pos.
    left_idx  = np.where(x_pos < left_range)[0]
    right_idx = np.where(x_pos >= left_range)[0]

    ordered_left  = left_idx[np.argsort(y_pos[left_idx])]
    ordered_right = right_idx[np.argsort(y_pos[right_idx])]
    ordered = list(ordered_left) + list(ordered_right)

    return [regions[i] for i in ordered]


# ── Stage 2a: stain binarisation (binarize_card.m) ───────────────────────────

def binarize(R_crop: np.ndarray, mask_crop: np.ndarray) -> np.ndarray:
    """binarize_card.m: stain mask = (R < threshold) & card_mask."""
    return (R_crop < STAIN_THRESHOLD) & mask_crop


# ── Stage 2b: component labelling (detect_droplets.m) ────────────────────────

def detect_components(bw: np.ndarray):
    """detect_droplets.m: bwconncomp(bw, 8) → regionprops.  No area filter applied."""
    lbl = label(bw, connectivity=2)   # connectivity=2 → 8-connectivity in 2-D
    return regionprops(lbl)


# ── Stage 2c: sub-segmentation (isolate_elements.m) ──────────────────────────

def isolate_elements(components, R_crop: np.ndarray,
                      bw_crop: np.ndarray | None = None):
    """isolate_elements.m: adaptive R-channel threshold sub-segmentation.

    For every component, crop the R channel to the component bounding box,
    compute an adaptive threshold, and run detect_droplets on the sub-binary
    image. If no sub-regions are found, keep the original component.

    MATLAB (isolate_elements.m):
        element = imcrop(R, components(idx).BoundingBox)
        threshold = 190 - (190 - min(element(:))) / 2
        [~, new_components, ~] = detect_droplets(element < threshold)
        if numel(new_components) == 0
            keep original component
        else
            keep all sub-components (no area filter)

    Two MATLAB-specific behaviours are reproduced here:
    1. imcrop extent: MATLAB imcrop(R, BoundingBox) returns one extra row and one extra
       column beyond the skimage bbox because BoundingBox endpoints are half-integer and
       MATLAB rounds them up.  The crop is therefore extended by +1 in each direction
       (clipped to image boundaries).
    2. Integer threshold: R is uint8 in MATLAB, so
       ``190 - (190 - min(element)) / 2`` uses integer division with round-half-to-even.
       Python's built-in round() replicates this; a floating-point division would include
       one extra pixel value for blobs whose (190 - min_val) is odd.

    Parameters
    ----------
    components : list of regionprops (or _OffsetRegion) objects
    R_crop     : red-channel card crop (float64 or uint8)
    bw_crop    : optional stain binary card crop (bool).  When provided, the
                 function also returns a reconstructed binary mask showing the
                 effective pixels after adaptive sub-segmentation (useful for
                 step-by-step comparison with MATLAB).  For components where
                 no sub-regions are found the original stain pixels are kept;
                 for components that are sub-segmented the adaptive-threshold
                 binary (element < threshold) is used instead.

    Returns
    -------
    all_components          : list of components (as before)
    isol_mask (only when bw_crop is provided) : bool ndarray, same shape as R_crop
    """
    all_components = []
    build_mask = bw_crop is not None
    if build_mask:
        isol_mask = np.zeros(R_crop.shape[:2], dtype=bool)

    h_crop, w_crop = R_crop.shape[:2]
    for comp in components:
        rb0, cb0, rb1, cb1 = comp.bbox   # skimage 0-indexed, exclusive end

        # MATLAB imcrop(R, BoundingBox) returns one extra row and one extra column
        # beyond the skimage bbox (round(ymin+height) and round(xmin+width) each
        # include the next pixel due to half-integer BoundingBox coordinates).
        # Replicate this by extending the crop by 1 in each direction, clipped to bounds.
        re1 = min(rb1 + 1, h_crop)
        ce1 = min(cb1 + 1, w_crop)
        element = R_crop[rb0:re1, cb0:ce1].astype(np.float64)

        # MATLAB R is uint8, so the threshold formula uses integer arithmetic:
        # threshold = 190 - (190 - min(element)) / 2  (uint8 division, round-half-to-even)
        # When (190 - min_val) is odd, floating-point gives a .5 threshold that
        # includes one extra pixel value compared to MATLAB's integer result.
        # Python's built-in round() also uses round-half-to-even, matching MATLAB.
        min_val   = int(element.min())
        threshold = ADAPTIVE_INITIAL - round((ADAPTIVE_INITIAL - min_val) / ADAPTIVE_FRACTION)

        sub_bw   = element < threshold
        sub_comp = detect_components(sub_bw)   # no area filter (mirrors MATLAB)

        if len(sub_comp) == 0:
            all_components.append(comp)
            if build_mask:
                # No sub-segmentation: keep original stain pixels for this component
                isol_mask[rb0:rb1, cb0:cb1] |= bw_crop[rb0:rb1, cb0:cb1]
        else:
            # Offset sub-component centroids to card-crop coordinates
            for sc in sub_comp:
                all_components.append(_OffsetRegion(sc, rb0, cb0))
            if build_mask:
                # Sub-segmented: use adaptive-threshold binary (may differ from stainmask)
                isol_mask[rb0:re1, cb0:ce1] |= sub_bw

    if build_mask:
        return all_components, isol_mask
    return all_components


class _OffsetRegion:
    """Thin wrapper that shifts centroid & bbox of a sub-region back to crop coords."""

    def __init__(self, region, row_offset: int, col_offset: int):
        self.area     = region.area
        r0, c0, r1, c1 = region.bbox
        self._centroid = (
            region.centroid[0] + row_offset,
            region.centroid[1] + col_offset,
        )
        self.bbox = (r0 + row_offset, c0 + col_offset, r1 + row_offset, c1 + col_offset)

    @property
    def centroid(self):
        return self._centroid


# ── Stage 3: statistics (droplet_statistics.m) ───────────────────────────────

def droplet_statistics(components, dpi: int, spread_factor: bool = False) -> dict:
    """droplet_statistics.m: area → diameter conversion, NMD/VMD/CH statistics.

    MATLAB std() uses N-1 normalisation (ddof=1 in NumPy).
    """
    pixel_um = 25_400.0 / dpi   # microns per pixel (MATLAB: 25.4/(dpi/1000))

    # Build area histogram: area_hist[k-1] = # droplets with area == k px
    if not components:
        return {}

    areas    = np.array([c.area for c in components], dtype=int)
    max_area = int(areas.max())
    hist     = np.zeros(max_area, dtype=int)
    for a in areas:
        hist[a - 1] += 1   # 0-indexed; hist[k-1] = count for area k

    # Diameters (with or without spread factor)
    idx_1based = np.arange(1, max_area + 1, dtype=np.float64)
    stain_diam = pixel_um * np.sqrt(4.0 * idx_1based / np.pi)
    if spread_factor:
        diam = 0.53549306 * stain_diam - 0.000084839 * stain_diam**2
    else:
        diam = stain_diam

    # Expand to per-droplet list (MATLAB: full_diametri)
    full_d = np.repeat(diam, hist)
    if len(full_d) == 0:
        return {}

    mean_d = float(np.mean(full_d))
    std_d  = float(np.std(full_d, ddof=1))   # N-1 (MATLAB std default)

    n = len(diam)

    def _pct(cumulative, frac):
        # MATLAB: diametri(1 + sum(cumsum <= frac*total)) is 1-indexed.
        # Python 0-indexed equivalent: diam[sum(cumulative <= frac*total)]
        idx = int(np.sum(cumulative <= frac * cumulative[-1]))
        return diam[min(idx, n - 1)]

    # NMD family (cumulative by count)
    cum_n  = np.cumsum(hist)
    nmd    = _pct(cum_n, 0.50)
    nmd10  = _pct(cum_n, 0.10)
    nmd90  = _pct(cum_n, 0.90)

    # VMD family (cumulative by volume)
    vol   = hist * ((4.0 / 3.0) * np.pi * (diam / 2.0)**3)
    cum_v = np.cumsum(vol)
    vmd   = _pct(cum_v, 0.50)
    vmd10 = _pct(cum_v, 0.10)
    vmd90 = _pct(cum_v, 0.90)

    ch = vmd / nmd if nmd > 0 else float("nan")

    # Cumulative count percentage (for the cumulative distribution plot)
    csn = 100.0 * cum_n.astype(float) / float(cum_n[-1])

    # Per-bin droplet histogram (for the grouped bar chart)
    n_bins = len(DIAMETER_LIMITS) + 1
    bin_idx = np.digitize(diam, DIAMETER_LIMITS)           # values 0 … n_bins-1
    droplet_histogram = np.bincount(
        bin_idx, weights=hist.astype(float), minlength=n_bins
    ).astype(int)

    return dict(
        nmd=nmd, nmd10=nmd10, nmd90=nmd90,
        vmd=vmd, vmd10=vmd10, vmd90=vmd90,
        ch=ch, mean=mean_d, std=std_d,
        diam=diam, hist=hist, csn=csn,
        droplet_histogram=droplet_histogram,
    )


def _print_results(card_labels, stats_raw_list, stats_sf_list, coverage_list, n_final_list):
    """Print MATLAB-style results table: metrics as rows, cards as columns.

    Mirrors the console output format of compute_droplet_coverage.m.
    """
    label_w = 22   # width of the row-label column
    val_w   = 9    # width of each value column

    header = f"{'Card analysis:':<{label_w}}" + ''.join(f"{lbl:>{val_w}}" for lbl in card_labels)
    sep    = '-' * len(header)

    def _row(label, values, fmt='.2f'):
        cells = ''.join(
            f"{(v if isinstance(v, str) else f'{v:{fmt}}'):>{val_w}}"
            for v in values
        )
        print(f"{label:<{label_w}}{cells}")

    def _get(stats_list, key):
        return [s.get(key, float('nan')) if s else float('nan') for s in stats_list]

    print()
    print(header)
    print(sep)
    print('Without Spread Factor')
    _row('NMD',           _get(stats_raw_list, 'nmd'))
    _row('NMD10',         _get(stats_raw_list, 'nmd10'))
    _row('NMD90',         _get(stats_raw_list, 'nmd90'))
    _row('VMD',           _get(stats_raw_list, 'vmd'))
    _row('VMD10',         _get(stats_raw_list, 'vmd10'))
    _row('VMD90',         _get(stats_raw_list, 'vmd90'))
    _row('CH',            _get(stats_raw_list, 'ch'))
    _row('Mean diameter', _get(stats_raw_list, 'mean'))
    _row('Diameter STD',  _get(stats_raw_list, 'std'))
    print(sep)
    print('With Spread Factor')
    _row('NMD',           _get(stats_sf_list, 'nmd'))
    _row('NMD10',         _get(stats_sf_list, 'nmd10'))
    _row('NMD90',         _get(stats_sf_list, 'nmd90'))
    _row('VMD',           _get(stats_sf_list, 'vmd'))
    _row('VMD10',         _get(stats_sf_list, 'vmd10'))
    _row('VMD90',         _get(stats_sf_list, 'vmd90'))
    _row('CH',            _get(stats_sf_list, 'ch'))
    _row('Mean diameter', _get(stats_sf_list, 'mean'))
    _row('Diameter STD',  _get(stats_sf_list, 'std'))
    print(sep)
    _row('Covered area',      [f'{v:.2f}%'  for v in coverage_list])
    _row('Droplets per card', [f'{n:>d}'    for n in n_final_list])
    drops_cm2 = [n / CARD_AREA_CM2 for n in n_final_list]
    _row('Droplets per cm2',  drops_cm2)
    print(sep)
    print()


# ── Output helpers (analyze_card_image.m) ────────────────────────────────────

def _image_name_from_path(path: str) -> str:
    """Derive plot/CSV filename stem — mirrors analyze_card_image.m naming.

    MATLAB splits on '.' then '_', joins parts with '-':
        '20250718_A_UAS_56.tif'  →  '20250718-A-UAS-56'
    """
    return Path(path).stem.replace('_', '-')


def _save_segmentation_plot(mask, cards, card_labels, image_name, results_folder):
    """Save card segmentation overlay — mirrors create_card_mask.m figure."""
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.imshow(mask, cmap='gray', interpolation='nearest')
    for card_region, lbl in zip(cards, card_labels):
        cy, cx = card_region.centroid   # regionprops: (row, col) → matplotlib (x=col, y=row)
        ax.text(cx, cy, lbl, color='red', fontsize=14,
                ha='center', va='center', fontweight='bold')
    ax.set_title(image_name)
    ax.axis('off')
    fig.savefig(Path(results_folder) / f'cards-{image_name}.png',
                dpi=150, bbox_inches='tight')
    plt.close(fig)


def _save_binary_mask(full_binary, image_name, results_folder):
    """Save binary stain mask — mirrors imwrite(full_binary, ...) in analyze_card_image.m."""
    plt.imsave(str(Path(results_folder) / f'mask-cards-{image_name}.png'),
               full_binary, cmap='gray')


def _save_histogram_plot(droplet_histograms, card_labels, image_name, results_folder):
    """Grouped bar histogram — mirrors hist-<name>.png in analyze_card_image.m."""
    n_cards = len(card_labels)
    n_bins  = len(DIAMETER_LIMITS) + 1

    categories = [f'<{DIAMETER_LIMITS[0]}']
    for i in range(len(DIAMETER_LIMITS) - 1):
        categories.append(f'{DIAMETER_LIMITS[i]}–{DIAMETER_LIMITS[i + 1]}')
    categories.append(f'>{DIAMETER_LIMITS[-1]}')

    x     = np.arange(n_bins)
    width = 0.8 / n_cards

    fig, ax = plt.subplots(figsize=(8, 6))
    for ci, (hist, lbl) in enumerate(zip(droplet_histograms, card_labels)):
        offset = (ci - n_cards / 2 + 0.5) * width
        ax.bar(x + offset, hist, width, label=lbl,
               color=COLORS[ci % len(COLORS)],
               edgecolor='black', linewidth=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(categories, rotation=45, ha='right')
    ax.set_xlabel('Diameter (µm)')
    ax.set_ylabel('Droplets')
    ax.set_title('Droplets Number')
    ax.legend(loc='upper right', bbox_to_anchor=(1.18, 1))
    fig.savefig(Path(results_folder) / f'hist-{image_name}.png',
                dpi=300, bbox_inches='tight')
    plt.close(fig)


def _save_cumulative_plot(stats_list_raw, card_labels, image_name, results_folder):
    """Cumulative diameter distribution — mirrors cumulative-<name>.png in analyze_card_image.m."""
    max_diam = 1200
    fig, ax  = plt.subplots(figsize=(8, 6))

    for stats, lbl, color in zip(stats_list_raw, card_labels, COLORS):
        if not stats:
            continue
        diam = stats['diam']
        csn  = stats['csn']
        keep = diam <= max_diam
        ax.plot(np.concatenate([[0.0], diam[keep]]),
                np.concatenate([[0.0], csn[keep]]),
                color=color, linewidth=1, label=lbl)

    ax.set_xlabel('Diameter (µm)')
    ax.set_ylabel('Droplets (%)')
    ax.set_title('Cumulative Sum', fontweight='bold', pad=10)
    ax.set_xlim(0, max_diam)
    ax.set_ylim(0, 100)
    ax.set_xticks(range(0, max_diam + 1, 100))
    ax.tick_params(axis='x', rotation=45)
    ax.legend(loc='lower right')

    # MATLAB box on — ensure all four spines are visible
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.8)

    fig.savefig(Path(results_folder) / f'cumulative-{image_name}.png',
                dpi=300, bbox_inches='tight', pad_inches=0.15)
    plt.close(fig)


def _save_csv(rows, image_name, results_folder):
    """Write per-card statistics CSV — mirrors scan_cards.m output format."""
    out_path = Path(results_folder) / f'{image_name}.csv'
    with open(out_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=CSV_HEADER)
        writer.writeheader()
        writer.writerows(rows)
    print(f"  CSV  -> {out_path}")


# ── Main pipeline ─────────────────────────────────────────────────────────────

def run(image_path: str, dpi: int = 600, results_folder: str | None = None,
        save_masks: bool = False, save_intermediate_masks: bool = False) -> None:
    """Full MATLAB-equivalent pipeline for one scan image.

    Replicates compute_droplet_coverage.m + analyze_card_image.m.
    Results are always printed to stdout.  When results_folder is given,
    four PNG plots and a CSV file are also saved there.
    When save_masks is True, per-card binary PNG masks are also saved
    (cardmask-WSP-N-*.png and stainmask-WSP-N-*.png) for direct comparison
    with MATLAB output from export_masks.m.

    The --save-intermediate-masks flag adds two further masks per card for
    pinpointing numerical differences vs MATLAB at each pipeline stage:
        step1mask-WSP-N-*.png  after detect_droplets   (labeled > 0; should equal stainmask)
        step2mask-WSP-N-*.png  after isolate_elements  (reconstructed from adaptive sub-seg)

    Note: binary_fill_holes and remove_small_objects are used only for card-mask
    building (create_card_mask.m) and are NOT part of the stain analysis pipeline.
    droplet_statistics applies no area filter, so step2mask is the final binary.

    Parameters
    ----------
    image_path              : path to the scan image
    dpi                     : scanner resolution in dots per inch (default 600)
    results_folder          : if provided, saves cards-*.png, mask-cards-*.png,
                              hist-*.png, cumulative-*.png, and <image_name>.csv
    save_masks              : if True, saves per-card cardmask and stainmask PNGs
                              to results_folder (or to a 'masks' subfolder next to
                              the image when results_folder is not set)
    save_intermediate_masks : if True, saves step1mask and step2mask PNGs (same
                              folder as save_masks output) for pipeline diagnostics
    """
    print(f"\nReading {image_path}")
    img = tifffile.imread(image_path)

    if img.ndim == 2:
        img = np.stack([img] * 3, axis=2)
    if img.shape[2] > 3:
        img = img[:, :, :3]

    print(f"  dtype={img.dtype}  shape={img.shape}  "
          f"min={img.min()}  max={img.max()}")

    R_full = img[:, :, 0].astype(np.float64)
    G_full = img[:, :, 1].astype(np.float64)
    B_full = img[:, :, 2].astype(np.float64)

    # ── Stage 1: card mask (create_card_mask.m) ───────────────────────────────
    print("\nStage 1: building card mask ...")
    mask = build_card_mask(R_full, G_full, B_full)

    lbl     = label(mask, connectivity=2)
    regions = regionprops(lbl)
    print(f"  {len(regions)} card region(s) found")

    cards      = order_cards(regions, mask.shape)
    image_name = _image_name_from_path(image_path)

    # Resolve folder for per-card mask PNGs
    masks_folder: Path | None = None
    if save_masks or save_intermediate_masks:
        masks_folder = Path(results_folder) if results_folder else \
                       Path(image_path).parent / "masks"
        masks_folder.mkdir(parents=True, exist_ok=True)

    # ── Stage 2 & 3: per-card analysis ───────────────────────────────────────
    print(f"\nAnalysis ({len(cards)} cards):\n")
    print(f"{'Card':<10} {'Drops':>6} {'Init':>6} {'Cover%':>8} "
          f"{'Drops/cm2':>10} {'VMD_SF':>8} {'NMD_SF':>8} {'CH_SF':>6}")
    print("-" * 68)

    card_labels      = [f"WSP-{i + 1}" for i in range(len(cards))]
    full_binary      = np.zeros(img.shape[:2], dtype=bool)
    stats_raw_list   = []
    stats_sf_list    = []
    droplet_hist_sf  = []
    coverage_list    = []
    n_final_list     = []
    csv_rows         = []

    for card_idx, card_region in enumerate(cards):
        card_label     = card_labels[card_idx]
        card_mask_full = (lbl == card_region.label)

        # MATLAB uses logical indexing: sum(mask,2)>0, sum(mask,1)>0
        row_has_px = np.any(card_mask_full, axis=1)
        col_has_px = np.any(card_mask_full, axis=0)

        R_crop    = R_full[row_has_px, :][:, col_has_px]
        mask_crop = card_mask_full[row_has_px, :][:, col_has_px]

        # binarize_card.m: bw = (R < 120) & card_mask
        bw = binarize(R_crop, mask_crop)

        slug = card_label.replace(' ', '-')   # used by both mask-saving branches

        # Save per-card binary masks for direct comparison with MATLAB output
        if save_masks:
            plt.imsave(str(masks_folder / f"cardmask-{slug}-{image_name}.png"),
                       mask_crop.astype(np.uint8) * 255, cmap='gray')
            plt.imsave(str(masks_folder / f"stainmask-{slug}-{image_name}.png"),
                       bw.astype(np.uint8) * 255, cmap='gray')

        # Map bw back to full-image coordinates for the stain mask PNG
        rows = np.where(row_has_px)[0]
        cols = np.where(col_has_px)[0]
        full_binary[np.ix_(rows, cols)] |= bw

        stain_px = int(bw.sum())
        mask_px  = int(mask_crop.sum())
        coverage = 100.0 * stain_px / mask_px if mask_px > 0 else 0.0

        # detect_droplets.m: bwconncomp(bw, 8) — no area filter
        init_comps = detect_components(bw)

        # isolate_elements.m: adaptive sub-segmentation — no area filter
        # When --save-intermediate-masks is set, also reconstruct the per-step
        # binary masks so MATLAB and Python output can be compared step by step:
        #   step1mask  after detect_droplets  (labeled > 0 — should equal stainmask)
        #   step2mask  after isolate_elements (adaptive sub-threshold reconstruction)
        if save_intermediate_masks:
            # step1mask: labeled > 0 reconstructed from detect_droplets output
            step1_lbl = label(bw, connectivity=2)
            plt.imsave(str(masks_folder / f"step1mask-{slug}-{image_name}.png"),
                       (step1_lbl > 0).astype(np.uint8) * 255, cmap='gray')
            # step2mask: reconstructed binary from isolate_elements adaptive sub-seg
            final_comps, isol_mask = isolate_elements(init_comps, R_crop, bw)
            plt.imsave(str(masks_folder / f"step2mask-{slug}-{image_name}.png"),
                       isol_mask.astype(np.uint8) * 255, cmap='gray')
        else:
            final_comps = isolate_elements(init_comps, R_crop)
        n_final     = len(final_comps)

        # droplet_statistics.m (raw and with spread-factor correction)
        stats_raw = droplet_statistics(final_comps, dpi, spread_factor=False)
        stats_sf  = droplet_statistics(final_comps, dpi, spread_factor=True)

        stats_raw_list.append(stats_raw)
        stats_sf_list.append(stats_sf)
        coverage_list.append(coverage)
        n_final_list.append(n_final)
        droplet_hist_sf.append(
            stats_sf.get('droplet_histogram', np.zeros(len(DIAMETER_LIMITS) + 1, dtype=int))
        )

        def _g(d, k):
            return round(float(d.get(k, float('nan'))), 4)

        csv_rows.append({
            'filename': Path(image_path).name,  'card': card_label,
            'NMD':   _g(stats_raw, 'nmd'),   'NMD10': _g(stats_raw, 'nmd10'), 'NMD90': _g(stats_raw, 'nmd90'),
            'VMD':   _g(stats_raw, 'vmd'),   'VMD10': _g(stats_raw, 'vmd10'), 'VMD90': _g(stats_raw, 'vmd90'),
            'CH':    _g(stats_raw, 'ch'),
            'mean_diameter': _g(stats_raw, 'mean'), 'std_diameter': _g(stats_raw, 'std'),
            'covered_area_pct': round(coverage, 4),
            'droplet_count': n_final,
            'droplets_per_cm2': round(n_final / CARD_AREA_CM2, 4),
            'NMD_SF':   _g(stats_sf, 'nmd'),   'NMD10_SF': _g(stats_sf, 'nmd10'), 'NMD90_SF': _g(stats_sf, 'nmd90'),
            'VMD_SF':   _g(stats_sf, 'vmd'),   'VMD10_SF': _g(stats_sf, 'vmd10'), 'VMD90_SF': _g(stats_sf, 'vmd90'),
            'CH_SF':    _g(stats_sf, 'ch'),
            'mean_diameter_SF': _g(stats_sf, 'mean'), 'std_diameter_SF': _g(stats_sf, 'std'),
        })

    # ── Print MATLAB-style table (all cards as columns) ───────────────────────
    _print_results(card_labels, stats_raw_list, stats_sf_list, coverage_list, n_final_list)

    # ── Save outputs (analyze_card_image.m) ───────────────────────────────────
    if results_folder:
        Path(results_folder).mkdir(parents=True, exist_ok=True)
        print(f"Saving results to {results_folder} ...")
        _save_segmentation_plot(mask, cards, card_labels, image_name, results_folder)
        print(f"  PNG  -> cards-{image_name}.png")
        _save_binary_mask(full_binary, image_name, results_folder)
        print(f"  PNG  -> mask-cards-{image_name}.png")
        _save_histogram_plot(droplet_hist_sf, card_labels, image_name, results_folder)
        print(f"  PNG  -> hist-{image_name}.png")
        _save_cumulative_plot(stats_raw_list, card_labels, image_name, results_folder)
        print(f"  PNG  -> cumulative-{image_name}.png")
        _save_csv(csv_rows, image_name, results_folder)
    if save_masks or save_intermediate_masks:
        print(f"Per-card masks saved to {masks_folder}")


# ── CLI ───────────────────────────────────────────────────────────────────────

def main():
    p = argparse.ArgumentParser(
        description="MATLAB-equivalent WSP analysis (DPI and optional output folder).",
    )
    p.add_argument("image", help="Path to the scan image (TIF/PNG/JPG …)")
    p.add_argument("--dpi", type=int, default=600, help="Scanner DPI (default: 600)")
    p.add_argument("--results-folder", default=None,
                   help="Save CSV and PNG plots to this folder (created if absent)")
    p.add_argument("--save-masks", action="store_true",
                   help="Save per-card binary cardmask and stainmask PNGs for "
                        "comparison with MATLAB export_masks.m output")
    p.add_argument("--save-intermediate-masks", action="store_true",
                   help="Save step1mask (after detect_droplets) and step2mask "
                        "(after isolate_elements) per card for pipeline diagnostics. "
                        "These reveal which stage causes MATLAB/Python differences.")
    args = p.parse_args()

    if not Path(args.image).exists():
        sys.exit(f"Error: file not found: {args.image}")

    run(args.image, dpi=args.dpi, results_folder=args.results_folder,
        save_masks=args.save_masks,
        save_intermediate_masks=args.save_intermediate_masks)


if __name__ == "__main__":
    main()
