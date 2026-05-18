# WSP Scan — Water Sensitive Paper Analysis

Automated analysis pipeline for **Water Sensitive Paper (WSP)** cards used to
characterise spray deposits from drone or field sprayers.

A flatbed scanner image containing multiple WSP cards is processed to extract,
for each card:

- droplet count and spatial density (droplets/cm²)
- covered area percentage
- numeric and volume median diameters (NMD, VMD) with 10th/90th percentile variants
- coefficient of homogeneity (CH)
- mean and standard deviation of droplet diameter

Results are reported both as raw stain diameters and corrected with a spread factor.

---

## Repository Layout

| Folder | Description |
|---|---|
| [`matlab/`](matlab/README.md) | Original MATLAB implementation — reference pipeline |
| [`matlab2python/`](matlab2python/README.md) | Direct line-by-line Python translation of the MATLAB pipeline |
| `python/` | Enhanced Python version *(WORK IN PROGRESS)* |

---

## How it Works

1. **Card segmentation** — WSP cards are located in the scan using an RGB colour
   threshold, morphological hole-filling, and disk closing.
2. **Stain binarisation** — the red channel of each card is thresholded to isolate
   dark stain marks from the yellow background.
3. **Droplet detection** — 8-connected components are extracted from the binary
   stain mask.
4. **Sub-segmentation** — overlapping droplet blobs are split using a per-blob
   adaptive threshold on the red channel.
5. **Statistics** — component area is converted to equivalent circular diameter and
   summary statistics are computed from the resulting diameter histogram.

---

## Where to Start

- To run the analysis today, use the **[`matlab2python/`](matlab2python/README.md)**
  version — it is fully functional, numerically equivalent to the original MATLAB
  code, and runs on Windows, macOS, and Linux.
- To read or modify the reference algorithm, see **[`matlab/`](matlab/README.md)**.
