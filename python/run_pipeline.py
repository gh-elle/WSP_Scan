"""run_pipeline.py — WSP droplet analysis pipeline entry point.

Accepts a single image file or a directory of image files and runs the full
three-stage pipeline on each image:

  Stage 1 — card_masks.py      : segment WSP cards → overlay + stages images
  Stage 2 — detect_droplets.py : detect droplets   → per-card detection images
  Stage 3 — droplet_statistics : compute statistics → CSV

Output files for each image are written to ``<output_dir>/<image_stem>/``.
If processing a directory, all per-image CSV files are also merged into a
single ``<output_dir>/statistics_all.csv``.

Usage
-----
    # single image
    python run_pipeline.py card_example.tif

    # directory of images
    python run_pipeline.py scans/

    # with overrides
    python run_pipeline.py scans/ --output results/ --config my_config.yaml \\
        --colour-space lab --yellow-b 18 --blue-b -10 --droplet-b 3 \\
        --border-erosion 20 --overlay-alpha 0.4

    # MATLAB-compatible RGB mode (colour space + sub-segmentation)
    python run_pipeline.py scans/ --colour-space rgb

    # Results identical to MATLAB (RGB, no border erosion, adaptive threshold, MATLAB closing SE):
    python run_pipeline.py scans/ --colour-space rgb --border-erosion 0 \\
        --subseg-method adaptive_threshold --closing-se matlab_disk
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd
from skimage.color import rgb2lab

from card_masks import (
    Config, load_config, load_image,
    _build_card_mask, _cleanup, _get_regions, _order_cards,
    save_overlay, save_mask_stages, save_lab_bands, save_rgb_channels,
)
from detect_droplets import detect_droplets, save_card_detections
from droplet_statistics import compute_card_stats, compute_plot_data, save_scan_plots, stats_to_dataframe, print_stats_table


_IMAGE_SUFFIXES = {".tif", ".tiff", ".png", ".jpg", ".jpeg", ".bmp"}


def process_image(
    image_path: Path,
    out_dir: Path,
    cfg: Config,
    card_names: list[str] | None,
) -> pd.DataFrame:
    """Run the full pipeline on a single scan image.

    Parameters
    ----------
    image_path:
        Path to the scan image file.
    out_dir:
        Directory where all output files for this image will be written.
    cfg:
        Pipeline configuration (already loaded and overrides applied).
    card_names:
        Optional custom WSP labels.

    Returns
    -------
    pd.DataFrame
        Statistics table for this image (one row per card).
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = image_path.stem

    print(f"\n{'='*60}")
    print(f"Processing: {image_path.name}")
    print(f"  Output  : {out_dir}")

    # ── Load ──────────────────────────────────────────────────────────────────
    rgb = load_image(image_path)
    print(f"  Size    : {rgb.shape[1]}×{rgb.shape[0]} px")
    print(f"  Mode    : {cfg.colour_space.upper()} colour space")

    # Convert to working image once; in LAB mode compute LAB, in RGB mode use RGB directly
    if cfg.colour_space == "lab":
        img = rgb2lab(rgb)
    else:
        img = rgb

    # ── Stage 1: card segmentation ────────────────────────────────────────────
    print("  Stage 1 : segmenting cards …")
    raw_mask   = _build_card_mask(img, cfg)
    clean_mask = _cleanup(raw_mask, cfg)
    regions    = _get_regions(clean_mask)
    cards      = _order_cards(regions, clean_mask, card_names, cfg)
    print(f"            {len(cards)} cards found")

    save_overlay(rgb, cards, out_dir / f"overlay_{stem}.png", cfg)
    save_mask_stages(rgb, raw_mask, clean_mask, cards,
                     out_dir / f"stages_{stem}.png", cfg)
    if cfg.colour_space == "lab":
        save_lab_bands(img, out_dir / f"lab_bands_{stem}.png", cfg)
    else:
        save_rgb_channels(rgb, out_dir / f"rgb_channels_{stem}.png", cfg)

    # ── Stage 2: droplet detection ────────────────────────────────────────────
    print("  Stage 2 : detecting droplets …")
    droplets_per_card = []
    stain_px_per_card = []
    for card in cards:
        drops, stain_px = detect_droplets(img, card, cfg)
        droplets_per_card.append(drops)
        stain_px_per_card.append(stain_px)
        print(f"            {card.label}: {len(drops)} droplets")

    save_card_detections(rgb, img, cards, droplets_per_card, cfg,
                         out_dir / f"droplets_{stem}.png")

    # ── Stage 3: statistics ───────────────────────────────────────────────────
    print("  Stage 3 : computing statistics …")
    stats_list    = []
    plot_data_list = []
    for card, drops, stain_px in zip(cards, droplets_per_card, stain_px_per_card):
        stats = compute_card_stats(
            card_name    = image_path.name,
            wsp_id       = card.label,
            droplets     = drops,
            mask_area_px = int(card.mask.sum()),
            stain_area_px= stain_px,
            dpi          = cfg.dpi,
        )
        stats_list.append(stats)
        plot_data_list.append(compute_plot_data(drops, cfg.dpi, cfg))
        print(f"            {card.label}: coverage {stats.covered_area_pct:.1f}%  "
              f"VMD {stats.vmd:.0f} µm (raw) / {stats.vmd_sf:.0f} µm (SF)")

    save_scan_plots(plot_data_list, stats_list,
                    out_dir / f"plots_{stem}", cfg)

    df = stats_to_dataframe(stats_list)
    csv_path = out_dir / f"statistics_{stem}.csv"
    df.to_csv(csv_path, index=False, float_format="%.2f")
    print(f"            CSV -> {csv_path}")

    print_stats_table(stats_list)

    return df


def main():
    p = argparse.ArgumentParser(
        description="WSP droplet analysis pipeline — process one image or a directory.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "input",
        help="Image file (TIF/PNG/JPG/…) or directory containing image files.",
    )
    p.add_argument(
        "--output", "-o", default=None,
        help="Output directory.  Default: <input_dir>/results  (or <input_file_dir>/results).",
    )
    p.add_argument(
        "--config", "-c", default=None, metavar="YAML",
        help="Config YAML file (default: config.yaml next to this script).",
    )
    p.add_argument(
        "--card-names", nargs="+", default=None, metavar="NAME",
        help="Custom WSP labels, applied to every image (e.g. A B C D E F).",
    )

    p.add_argument(
        "--colour-space", choices=["lab", "rgb"], default=None,
        help="Colour space for detection: 'lab' (default) or 'rgb' (MATLAB-compatible).",
    )

    # LAB parameter overrides
    p.add_argument("--yellow-b",      type=float, default=None,
                   help="LAB b* minimum for yellow paper detection.")
    p.add_argument("--blue-b",        type=float, default=None,
                   help="LAB b* maximum for blue stain detection.")
    p.add_argument("--min-chroma",    type=float, default=None,
                   help="Minimum LAB chroma for blue stain detection.")
    p.add_argument("--droplet-b",     type=float, default=None,
                   help="LAB b* upper bound for droplet binarisation.")
    p.add_argument("--droplet-l-max", type=float, default=None,
                   help="LAB L* upper bound for droplet binarisation (excludes bright scratches).")

    # RGB parameter overrides
    p.add_argument("--yellow-rg",     type=float, default=None,
                   help="RGB offset: yellow if R > B+x AND G > B+x (RGB mode).")
    p.add_argument("--blue-rg",       type=float, default=None,
                   help="RGB offset: blue if B > R+x AND B > G+x (RGB mode).")
    p.add_argument("--droplet-r",     type=float, default=None,
                   help="RGB R-channel threshold for droplet binarisation (RGB mode).")

    # Shared overrides
    p.add_argument("--border-erosion", type=int,  default=None,
                   help="Card-edge erosion in pixels (excludes border artefacts).")
    p.add_argument("--overlay-alpha",  type=float, default=None,
                   help="Opacity of the red stain overlay in detection images (0–1).")

    # Sub-segmentation overrides
    p.add_argument("--subseg-method", choices=["watershed", "adaptive_threshold"],
                   default=None,
                   help="Sub-segmentation method: 'watershed' (default) or "
                        "'adaptive_threshold' (MATLAB-compatible, RGB mode only).")
    p.add_argument("--adaptive-thresh", type=float, default=None,
                   metavar="THRESH",
                   help="Initial R-channel threshold for adaptive sub-segmentation (RGB mode).")
    p.add_argument("--adaptive-frac", type=float, default=None,
                   metavar="FRAC",
                   help="Fraction divisor for adaptive threshold (must be > 1).")
    p.add_argument("--closing-se", choices=["disk", "matlab_disk"], default=None,
                   help="Closing structuring element: 'disk' (true Euclidean, default) or "
                        "'matlab_disk' (MATLAB-compatible 16-sided polygon, ~1 px wider).")

    args = p.parse_args()

    input_path = Path(args.input)
    if not input_path.exists():
        sys.exit(f"Error: path not found: {input_path}")

    # Collect image files
    if input_path.is_file():
        if input_path.suffix.lower() not in _IMAGE_SUFFIXES:
            sys.exit(f"Error: unsupported file type: {input_path.suffix}")
        image_files = [input_path]
        default_out = input_path.parent / "results"
    else:
        image_files = sorted(
            f for f in input_path.iterdir()
            if f.is_file() and f.suffix.lower() in _IMAGE_SUFFIXES
        )
        if not image_files:
            sys.exit(f"Error: no image files found in {input_path}")
        default_out = input_path / "results"

    out_root = Path(args.output) if args.output else default_out

    # Load config once; apply CLI overrides
    cfg = load_config(Path(args.config) if args.config else None)
    if args.colour_space   is not None: cfg.colour_space    = args.colour_space
    if args.yellow_b       is not None: cfg.yellow_b_min    = args.yellow_b
    if args.blue_b         is not None: cfg.blue_b_max      = args.blue_b
    if args.min_chroma     is not None: cfg.min_chroma      = args.min_chroma
    if args.droplet_b      is not None: cfg.droplet_b_max   = args.droplet_b
    if args.droplet_l_max  is not None: cfg.droplet_l_max   = args.droplet_l_max
    if args.yellow_rg      is not None: cfg.yellow_rg_over_b = args.yellow_rg
    if args.blue_rg        is not None: cfg.blue_b_over_rg  = args.blue_rg
    if args.droplet_r      is not None: cfg.droplet_r_max   = args.droplet_r
    if args.border_erosion is not None: cfg.border_erosion_px = args.border_erosion
    if args.overlay_alpha  is not None: cfg.overlay_alpha     = args.overlay_alpha
    if args.subseg_method  is not None: cfg.subseg_method     = args.subseg_method
    if args.adaptive_thresh is not None: cfg.adaptive_initial_threshold = args.adaptive_thresh
    if args.adaptive_frac  is not None: cfg.adaptive_fraction = args.adaptive_frac
    if args.closing_se     is not None: cfg.closing_se        = args.closing_se

    print(f"WSP pipeline — {len(image_files)} image(s) to process")

    all_dfs = []
    for image_path in image_files:
        # Each image gets its own sub-directory so outputs never collide
        out_dir = out_root / image_path.stem if len(image_files) > 1 else out_root
        df = process_image(image_path, out_dir, cfg, args.card_names)
        all_dfs.append(df)

    # Merge statistics across all images when processing a directory
    if len(all_dfs) > 1:
        merged     = pd.concat(all_dfs, ignore_index=True)
        merged_csv = out_root / "statistics_all.csv"
        out_root.mkdir(parents=True, exist_ok=True)
        merged.to_csv(merged_csv, index=False, float_format="%.2f")
        print(f"\nMerged statistics ({len(merged)} rows) → {merged_csv}")

    print("\nAll done.")


if __name__ == "__main__":
    main()
