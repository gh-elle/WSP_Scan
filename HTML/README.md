# WSP Scan — Browser Web App

A self-contained web application that runs the full WSP droplet analysis pipeline
entirely in the browser — no Python, no MATLAB, no installation required.

**→ [Open the app](https://gh-elle.github.io/WSP_Scan/)**

---

## Functionalities

- Drag-and-drop or click-to-select scan images (TIF, PNG, JPEG)
- Configurable DPI (300 / 600 / 1200 or custom value)
- Step-by-step progress log in the page while analysis runs
- Per-card results table with all statistics (NMD, VMD, NMD₁₀/₉₀, VMD₁₀/₉₀, CH,
  mean, std, coverage %, droplets, droplets/cm²) — raw stain and spread-factor corrected
- Visualisations: card overlay (greyscale + bounding boxes + card labels),
  stain mask (black on white), diameter histogram, cumulative distribution
- Export results as CSV, TXT, and PNG plots

### Export — saving files

| Browser | Behaviour |
|---|---|
| Chrome / Edge | "Select results folder" button lets you choose any local folder (e.g. `HTML/results/`) — all files are saved there in one click |
| Firefox | File System Access API not supported — each file is downloaded individually to the browser's downloads folder |

---

## Limits with Respect to the MATLAB / Python Code

| Aspect | Web app | MATLAB / Python |
|---|---|---|
| **Processing speed** | Slower for large TIF files (>50 MB); runs on a single thread in the browser | Fast, multi-core capable |
| **TIF / Chart support** | Libraries are self-hosted — works fully offline | No network needed |
| **Batch processing** | One file at a time | Scriptable over many files |
| **Folder output** | Chrome/Edge: folder picker; Firefox: individual downloads | Full control over output path |
| **Spread factor** | USDA-ARS polynomial (same as DepositScan software) | Same |
| **Numerical results** | Equivalent to MATLAB/Python within floating-point precision | Reference implementation |
