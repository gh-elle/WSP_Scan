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

Results are reported both as raw stain diameters and corrected with the USDA-ARS spread-factor polynomial (the same equation used by the DepositScan software).

---

## Repository Layout

| Folder | Description |
|---|---|
| [`matlab/`](matlab/README.md) | Original MATLAB implementation — reference pipeline |
| [`matlab2python/`](matlab2python/README.md) | Direct line-by-line Python translation of the MATLAB pipeline |
| [`HTML/`](HTML/README.md) | Browser-based web app — runs the full pipeline client-side, no installation needed · **[Open app](https://gh-elle.github.io/WSP_Scan/)** |
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
- To run the analysis **in your browser** without any installation, open
  **[the web app](https://gh-elle.github.io/WSP_Scan/)**
  — the full pipeline runs client-side in
  JavaScript, fully offline (Chrome or Edge recommended for folder export). See **[`HTML/README.md`](HTML/README.md)**
  for features and limits.
- To read or modify the reference algorithm, see **[`matlab/`](matlab/README.md)**.

---

## References

The MATLAB pipeline implemented in this repository was developed and validated in the following works:

1. Ghiani, L., Sassu, A., Piccirilli, D., Marcialis, G. L., & Gambella, F. (2019, September). *Development of a matlab code for the evaluation of spray distribution with water-sensitive paper.* In International Mid-Term Conference of the Italian Association of Agricultural Engineering (pp. 845–853). Cham: Springer International Publishing.

2. Sassu, A., Psiroukis, V., Bettucci, F., Ghiani, L., Fountas, S., & Gambella, F. (2024). *Unmanned aerial system plant protection products spraying performance evaluation on a vineyard.* Precision Agriculture, 25(4), 2082–2112.

