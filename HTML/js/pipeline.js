'use strict';
/**
 * pipeline.js — WSP scan analysis pipeline (browser port of wsp_matlab.py).
 *
 * Replicates the MATLAB WSP pipeline:
 *   create_card_mask.m → binarize_card.m → detect_droplets.m →
 *   isolate_elements.m → droplet_statistics.m
 *
 * Key MATLAB-matching behaviours:
 *   - Banker's rounding in isolate_elements adaptive threshold
 *   - imcrop +1 row/col extension in isolate_elements
 *   - N−1 (ddof=1) standard deviation
 *   - Exact MATLAB strel('disk',50) loaded from matlab_strel_disk50.png
 *
 * All stages run at full resolution, mirroring Python/MATLAB exactly.
 */

// ── Constants (matching wsp_matlab.py) ────────────────────────────────────────
const CARD_AREA_CM2     = 2.6 * 7.6;
const STAIN_THRESHOLD   = 120;
const ADAPTIVE_INITIAL  = 190;
const ADAPTIVE_FRACTION = 2;
const MIN_CARD_AREA     = 1000;
const DIAMETER_LIMITS   = [50, 100, 150, 200, 300, 400, 500, 600];
const COLORS_HEX        = [
    '#ff0000', '#0000ff', '#00ff00', '#00ffff', '#ff00ff',
    '#ff8000', '#8080ff', '#ff8080', '#ffff00', '#80ff80', '#4dff4d', '#808080',
];


// ── Structuring element ────────────────────────────────────────────────────────

/**
 * The exact MATLAB strel('disk',50) kernel embedded as a base64 PNG.
 * Using an inline data URL guarantees the SE loads correctly regardless of
 * whether the page is served from file:// or an HTTP server.
 */
const SE_PNG_B64 =
    'data:image/png;base64,' +
    'iVBORw0KGgoAAAANSUhEUgAAAMgAAADICAAAAACIM/FCAAAAB3RJTUUH6gUQFAMTV05TGgAAAX1J' +
    'REFUeJzt3MFVRDEMA8CF/nuGGtbWf2vEzN2OdE/yegEAAAAAAAAAAMB/9JVd9/Oxk7Pr3uiRPjq6' +
    '7a0e4bOTy97skT08uOvtHtHTc6sGPZLHxzaNegTPTy0a9sgFCO0Z94glyKxZ9EhFiGxZ9QhlSCxZ' +
    '9siECOxY94ik2K8I9EjEWG+I9Ajk2C4I9dgHWc7HeqyT7MaDPbZRVtPRHsssm+Fwj12YxWy8xyrN' +
    'fPSBHps448lHeizyTAcf6jEPNJx7rMc40WzswR7TSKOpR3sMM02GHu4xCzWYebzHKNV3PsVnKHKN' +
    'Itcoco0i1yhyjSLXKHKNItcoco0i1yhyjSLXKHKNItcoco0i1yhyjSLXKHKNItfUFPnPd1F6bgf1' +
    '3NfquUHXc6ex55Zpz73fnpvYPXfje14r9Lwf6XnR0/PGqufVW887xJ6XoT1vdXteT/e8Z+/5YaDn' +
    'z4eeXzh6/kXp+amm5++gnt+cev7X+uSPZwAAAAAAAAAAAMCf9gtQ2jmOtOZCfAAAAABJRU5ErkJg' +
    'gg==';

/**
 * Load the exact MATLAB strel('disk',50) kernel from the embedded PNG data URL.
 * Returns {se, seH, seW, cy, cx} where se is a Uint8Array of the bounding-box crop.
 */
function loadMatlabSE() {
    return new Promise((resolve, reject) => {
        const img = new Image();
        img.onload = () => {
            const cv  = document.createElement('canvas');
            cv.width  = img.width;
            cv.height = img.height;
            const ctx = cv.getContext('2d');
            ctx.drawImage(img, 0, 0);
            const px = ctx.getImageData(0, 0, img.width, img.height).data;

            let r0 = img.height, c0 = img.width, r1 = -1, c1 = -1;
            for (let r = 0; r < img.height; r++)
                for (let c = 0; c < img.width; c++)
                    if (px[(r * img.width + c) * 4] > 128) {
                        if (r < r0) r0 = r;  if (r > r1) r1 = r;
                        if (c < c0) c0 = c;  if (c > c1) c1 = c;
                    }
            if (r1 < 0) { reject(new Error('No white pixels found in SE image')); return; }

            const seH = r1 - r0 + 1, seW = c1 - c0 + 1;
            const se  = new Uint8Array(seH * seW);
            for (let r = r0; r <= r1; r++)
                for (let c = c0; c <= c1; c++)
                    if (px[(r * img.width + c) * 4] > 128)
                        se[(r - r0) * seW + (c - c0)] = 1;

            resolve({ se, seH, seW, cy: (seH - 1) >> 1, cx: (seW - 1) >> 1 });
        };
        img.onerror = () => reject(new Error('Failed to decode embedded SE PNG'));
        img.src = SE_PNG_B64;
    });
}

/** Downsample SE by factor using max-pooling (any SE pixel in a block → 1). */
function downsampleSE({ se, seH, seW }, factor) {
    const nH = Math.ceil(seH / factor), nW = Math.ceil(seW / factor);
    const ns = new Uint8Array(nH * nW);
    for (let r = 0; r < seH; r++)
        for (let c = 0; c < seW; c++)
            if (se[r * seW + c])
                ns[Math.floor(r / factor) * nW + Math.floor(c / factor)] = 1;
    return { se: ns, seH: nH, seW: nW, cy: (nH - 1) >> 1, cx: (nW - 1) >> 1 };
}

/**
 * Pre-compute horizontal row spans for the fast row-decomposition dilation.
 * Returns [{oy, ox_lo, ox_hi}] — one entry per non-empty SE row.
 * oy/ox_lo/ox_hi are offsets relative to the SE centre (cy, cx).
 */
function computeSESpans({ se, seH, seW, cy, cx }) {
    const spans = [];
    for (let dy = 0; dy < seH; dy++) {
        let lo = Infinity, hi = -Infinity;
        for (let dx = 0; dx < seW; dx++)
            if (se[dy * seW + dx]) {
                const ox = dx - cx;
                if (ox < lo) lo = ox;
                if (ox > hi) hi = ox;
            }
        if (lo <= hi) spans.push({ oy: dy - cy, ox_lo: lo, ox_hi: hi });
    }
    return spans;
}

/** Octagon SE approximation of MATLAB strel('disk', r, 4) — fallback only. */
function makeDiskSE(r) {
    const n = 4, size = 2 * r + 5, c = r + 2;
    const seH = size, seW = size;
    const se  = new Uint8Array(seH * seW);
    for (let row = 0; row < seH; row++)
        for (let col = 0; col < seW; col++) {
            const y = row - c, x = col - c;
            let inside = true;
            for (let k = 0; k < n; k++) {
                const angle = k * Math.PI / n;
                if (Math.abs(x * Math.cos(angle) + y * Math.sin(angle)) > r) { inside = false; break; }
            }
            if (inside) se[row * seW + col] = 1;
        }
    return { se, seH, seW, cy: c, cx: c };
}


// ── Morphological operations ──────────────────────────────────────────────────

/**
 * Binary dilation using row-decomposition with prefix sums.
 * For each SE row span {oy, ox_lo, ox_hi}, slides a range-OR window
 * over the corresponding source row.  O(|spans| × H × W).
 */
function dilate(src, spans, H, W) {
    const res    = new Uint8Array(H * W);
    const prefix = new Int32Array(W + 1);   // reused across all rows
    for (const { oy, ox_lo, ox_hi } of spans) {
        for (let rs = 0; rs < H; rs++) {
            const ro = rs + oy;
            if (ro < 0 || ro >= H) continue;

            // Build prefix sum for source row rs
            const sBase = rs * W;
            prefix[0] = 0;
            for (let c = 0; c < W; c++)
                prefix[c + 1] = prefix[c] + (src[sBase + c] ? 1 : 0);

            // For each output column, test if any src column in [c−ox_hi, c−ox_lo] is 1
            const oBase = ro * W;
            for (let c = 0; c < W; c++) {
                if (res[oBase + c]) continue;
                const lo = c - ox_hi > 0 ? c - ox_hi : 0;
                const hi = c - ox_lo + 1 < W ? c - ox_lo + 1 : W;
                if (lo < hi && prefix[hi] > prefix[lo]) res[oBase + c] = 1;
            }
        }
    }
    return res;
}

/** Binary erosion via complement-dilation (valid for symmetric SE). */
function erode(src, spans, H, W) {
    const inv = new Uint8Array(H * W);
    for (let i = 0; i < H * W; i++) inv[i] = src[i] ? 0 : 1;
    const d   = dilate(inv, spans, H, W);
    const res = new Uint8Array(H * W);
    for (let i = 0; i < H * W; i++) res[i] = d[i] ? 0 : 1;
    return res;
}

/** Morphological closing = dilate then erode. */
function morphClose(src, spans, H, W) {
    return erode(dilate(src, spans, H, W), spans, H, W);
}


// ── Binary image operations ────────────────────────────────────────────────────

/**
 * Fill holes: background pixels not 8-connected to the image border
 * are set to foreground.  Mirrors scipy.ndimage.binary_fill_holes.
 */
function binaryFillHoles(mask, H, W) {
    const vis = new Uint8Array(H * W);
    const qR  = new Int32Array(H * W);
    const qC  = new Int32Array(H * W);
    let head  = 0, tail = 0;

    const seed = (r, c) => {
        const i = r * W + c;
        if (!mask[i] && !vis[i]) { vis[i] = 1; qR[tail] = r; qC[tail++] = c; }
    };
    for (let r = 0; r < H; r++) { seed(r, 0); seed(r, W - 1); }
    for (let c = 1; c < W - 1; c++) { seed(0, c); seed(H - 1, c); }

    const DR = [-1, -1, -1, 0, 0, 1, 1, 1];
    const DC = [-1,  0,  1,-1, 1,-1, 0, 1];
    while (head < tail) {
        const r = qR[head], c = qC[head++];
        for (let d = 0; d < 8; d++) {
            const nr = r + DR[d], nc = c + DC[d];
            if (nr < 0 || nr >= H || nc < 0 || nc >= W) continue;
            const ni = nr * W + nc;
            if (!mask[ni] && !vis[ni]) { vis[ni] = 1; qR[tail] = nr; qC[tail++] = nc; }
        }
    }

    const res = new Uint8Array(H * W);
    for (let i = 0; i < H * W; i++) res[i] = (mask[i] || !vis[i]) ? 1 : 0;
    return res;
}

/**
 * 8-connectivity connected component labeling (BFS).
 * Returns {labels: Int32Array, regions: Array<{label, area, centroid, bbox}>}.
 * bbox = [r0, c0, r1+1, c1+1] (exclusive end, skimage convention).
 */
function labelCC(mask, H, W) {
    const labels = new Int32Array(H * W);
    const qR     = new Int32Array(H * W);
    const qC     = new Int32Array(H * W);
    const DR     = [-1, -1, -1, 0, 0, 1, 1, 1];
    const DC     = [-1,  0,  1,-1, 1,-1, 0, 1];
    const regions = [];
    let nextLabel = 0;

    for (let r0 = 0; r0 < H; r0++) {
        for (let c0 = 0; c0 < W; c0++) {
            const i0 = r0 * W + c0;
            if (!mask[i0] || labels[i0]) continue;

            const lbl = ++nextLabel;
            labels[i0] = lbl;
            let head = 0, tail = 0;
            qR[tail] = r0; qC[tail++] = c0;

            let area = 0, sumR = 0, sumC = 0;
            let bR0 = r0, bC0 = c0, bR1 = r0, bC1 = c0;
            while (head < tail) {
                const r = qR[head], c = qC[head++];
                area++; sumR += r; sumC += c;
                if (r < bR0) bR0 = r;  if (r > bR1) bR1 = r;
                if (c < bC0) bC0 = c;  if (c > bC1) bC1 = c;
                for (let d = 0; d < 8; d++) {
                    const nr = r + DR[d], nc = c + DC[d];
                    if (nr < 0 || nr >= H || nc < 0 || nc >= W) continue;
                    const ni = nr * W + nc;
                    if (mask[ni] && !labels[ni]) {
                        labels[ni] = lbl;
                        qR[tail] = nr; qC[tail++] = nc;
                    }
                }
            }
            regions.push({
                label:    lbl,
                area,
                centroid: [sumR / area, sumC / area],
                bbox:     [bR0, bC0, bR1 + 1, bC1 + 1],
            });
        }
    }
    return { labels, regions };
}

/** Remove connected components with area < minSize (8-connectivity). */
function removeSmallObjects(mask, H, W, minSize) {
    const { labels, regions } = labelCC(mask, H, W);
    const keep = new Uint8Array(regions.length + 1);
    for (const r of regions) if (r.area >= minSize) keep[r.label] = 1;
    const res = new Uint8Array(H * W);
    for (let i = 0; i < H * W; i++) if (labels[i] && keep[labels[i]]) res[i] = 1;
    return res;
}


// ── Stage 1: card mask ────────────────────────────────────────────────────────

/**
 * Replicate create_card_mask.m at full resolution — mirrors Python build_card_mask():
 *   colour threshold → fill holes → remove small objects → morph close.
 */
function buildCardMask(R, G, B, H, W, spans) {
    const mask = new Uint8Array(H * W);
    for (let i = 0; i < H * W; i++) {
        const r = R[i], g = G[i], b = B[i];
        mask[i] = ((r > b + 50 && g > b + 50) || (b > r + 15 && b > g + 15)) ? 1 : 0;
    }
    const filled = binaryFillHoles(mask, H, W);
    const opened = removeSmallObjects(filled, H, W, MIN_CARD_AREA);
    return morphClose(opened, spans, H, W);
}

/**
 * Replicate card ordering from create_card_mask.m:
 *   left column (col < W/3) cards sorted by row, then right-column cards by row.
 */
function orderCards(regions, W_ds) {
    const leftRange = Math.floor(W_ds / 3);
    const xPos      = regions.map(r => r.bbox[1]);   // col_min
    const yPos      = regions.map(r => r.bbox[0]);   // row_min

    const leftIdx  = regions.map((_, i) => i).filter(i => xPos[i] <  leftRange);
    const rightIdx = regions.map((_, i) => i).filter(i => xPos[i] >= leftRange);
    leftIdx.sort( (a, b) => yPos[a] - yPos[b]);
    rightIdx.sort((a, b) => yPos[a] - yPos[b]);

    return [...leftIdx, ...rightIdx].map(i => regions[i]);
}


// ── Stage 2: binarise & sub-segment ──────────────────────────────────────────

/** binarize_card.m: stain mask = (R < 120) & card_mask. */
function binarize(R_crop, mask_crop, size) {
    const bw = new Uint8Array(size);
    for (let i = 0; i < size; i++)
        bw[i] = (R_crop[i] < STAIN_THRESHOLD && mask_crop[i]) ? 1 : 0;
    return bw;
}

/**
 * Banker's rounding (round-half-to-even) — matches Python's round() and
 * MATLAB's uint8 integer division used in isolate_elements.m.
 */
function bankersRound(x) {
    const f    = Math.floor(x);
    const diff = x - f;
    if (Math.abs(diff - 0.5) < 1e-10) return (f % 2 === 0) ? f : f + 1;
    return Math.round(x);
}

/**
 * Replicate isolate_elements.m: adaptive R-channel threshold sub-segmentation.
 *
 * Two MATLAB-specific behaviours (see wsp_matlab.py for full explanation):
 *   1. imcrop +1 extension: bbox is extended by one extra row and column.
 *   2. Integer threshold:   bankersRound replicates MATLAB uint8 division.
 *
 * @param {Array}        comps    — array of {area, centroid, bbox} in crop coords
 * @param {Float32Array} R_crop   — red channel of the card crop
 * @param {number}       H_crop   — crop height
 * @param {number}       W_crop   — crop width
 * @returns {Array}               — final component list (same schema as input)
 */
function isolateElements(comps, R_crop, H_crop, W_crop) {
    const all = [];
    for (const comp of comps) {
        const [rb0, cb0, rb1, cb1] = comp.bbox;

        // MATLAB imcrop +1 extension (clipped to crop bounds)
        const re1 = rb1 + 1 < H_crop ? rb1 + 1 : H_crop;
        const ce1 = cb1 + 1 < W_crop ? cb1 + 1 : W_crop;
        const eH  = re1 - rb0, eW = ce1 - cb0;

        // Minimum value inside the sub-crop (element)
        let minVal = 255;
        for (let r = 0; r < eH; r++)
            for (let c = 0; c < eW; c++) {
                const v = R_crop[(rb0 + r) * W_crop + (cb0 + c)];
                if (v < minVal) minVal = v;
            }

        // Adaptive threshold with banker's rounding (matches MATLAB uint8 arithmetic)
        const threshold = ADAPTIVE_INITIAL - bankersRound((ADAPTIVE_INITIAL - minVal) / ADAPTIVE_FRACTION);

        // Sub-binarise
        const subBw = new Uint8Array(eH * eW);
        for (let r = 0; r < eH; r++)
            for (let c = 0; c < eW; c++)
                subBw[r * eW + c] = (R_crop[(rb0 + r) * W_crop + (cb0 + c)] < threshold) ? 1 : 0;

        const { regions: subComps } = labelCC(subBw, eH, eW);

        if (subComps.length === 0) {
            all.push(comp);
        } else {
            for (const sc of subComps) {
                all.push({
                    area:     sc.area,
                    centroid: [sc.centroid[0] + rb0, sc.centroid[1] + cb0],
                    bbox:     [sc.bbox[0] + rb0, sc.bbox[1] + cb0,
                               sc.bbox[2] + rb0, sc.bbox[3] + cb0],
                });
            }
        }
    }
    return all;
}


// ── Stage 3: statistics ───────────────────────────────────────────────────────

/** Return bin index matching np.digitize(d, DIAMETER_LIMITS). */
function digitizeBin(d) {
    let i = 0;
    while (i < DIAMETER_LIMITS.length && d >= DIAMETER_LIMITS[i]) i++;
    return i;
}

/**
 * Replicate droplet_statistics.m:
 *   area → stain diameter → (optionally) spread-factor diameter →
 *   NMD/VMD/CH/mean/std and histogram arrays.
 *
 * Uses N−1 (ddof=1) std to match MATLAB default.
 */
function dropletStatistics(comps, dpi, spreadFactor) {
    if (!comps.length) return null;
    const pixelUm = 25400.0 / dpi;
    const N       = comps.length;

    let maxArea = 0;
    for (const c of comps) if (c.area > maxArea) maxArea = c.area;

    const hist = new Int32Array(maxArea);
    for (const c of comps) hist[c.area - 1]++;

    // Diameter for each area bin (1-indexed areas)
    const diam = new Float64Array(maxArea);
    for (let k = 0; k < maxArea; k++) {
        const sd = pixelUm * Math.sqrt(4.0 * (k + 1) / Math.PI);
        diam[k] = spreadFactor ? 0.53549306 * sd - 0.000084839 * sd * sd : sd;
    }

    // Weighted mean and std (ddof=1)
    let meanD = 0;
    for (let k = 0; k < maxArea; k++) meanD += diam[k] * hist[k];
    meanD /= N;
    let varSum = 0;
    for (let k = 0; k < maxArea; k++) varSum += hist[k] * (diam[k] - meanD) ** 2;
    const stdD = N > 1 ? Math.sqrt(varSum / (N - 1)) : 0;

    // Cumulative count — used for NMD percentiles
    const cumN = new Float64Array(maxArea);
    cumN[0] = hist[0];
    for (let k = 1; k < maxArea; k++) cumN[k] = cumN[k - 1] + hist[k];

    // Cumulative volume — used for VMD percentiles
    const cumV = new Float64Array(maxArea);
    cumV[0] = hist[0] * (4 / 3) * Math.PI * (diam[0] / 2) ** 3;
    for (let k = 1; k < maxArea; k++)
        cumV[k] = cumV[k - 1] + hist[k] * (4 / 3) * Math.PI * (diam[k] / 2) ** 3;

    // Percentile helper — matches Python's _pct (counts elements <= frac*total)
    const pct = (cum, frac) => {
        const total = cum[maxArea - 1];
        let idx = 0;
        while (idx < maxArea && cum[idx] <= frac * total) idx++;
        return diam[Math.min(idx, maxArea - 1)];
    };

    const nmd   = pct(cumN, 0.50), nmd10 = pct(cumN, 0.10), nmd90 = pct(cumN, 0.90);
    const vmd   = pct(cumV, 0.50), vmd10 = pct(cumV, 0.10), vmd90 = pct(cumV, 0.90);
    const ch    = nmd > 0 ? vmd / nmd : NaN;

    // Cumulative % by count (for the distribution plot)
    const csn = new Float64Array(maxArea);
    for (let k = 0; k < maxArea; k++) csn[k] = 100.0 * cumN[k] / cumN[maxArea - 1];

    // Droplet count per diameter bin (for the grouped bar chart)
    const nBins          = DIAMETER_LIMITS.length + 1;
    const dropletHistogram = new Int32Array(nBins);
    for (let k = 0; k < maxArea; k++)
        dropletHistogram[digitizeBin(diam[k])] += hist[k];

    return { nmd, nmd10, nmd90, vmd, vmd10, vmd90, ch, mean: meanD, std: stdD,
             diam, hist, csn, dropletHistogram };
}


/**
 * Run the full WSP analysis pipeline on an RGBA image.
 *
 * @param {Uint8Array|Uint8ClampedArray} rgba  flat RGBA pixels, row-major (H×W×4)
 * @param {number}   H            image height in pixels
 * @param {number}   W            image width in pixels
 * @param {number}   dpi          scanner resolution (dots per inch)
 * @param {function} onProgress   called with a status string at each stage
 * @returns {Promise<object>}     analysis result object
 */
async function runPipeline(rgba, H, W, dpi, onProgress) {
    const yield_ = () => new Promise(r => setTimeout(r, 0));

    // ── 1. Load structuring element ───────────────────────────────────────────
    onProgress('1) Loading structuring element…');
    await yield_();
    let seData;
    try {
        seData = await loadMatlabSE();
        onProgress('   ✓ Structuring element loaded (exact MATLAB disk SE).');
    } catch (e) {
        onProgress(`   ⚠ ${e.message} — using octagon approximation.`);
        seData = makeDiskSE(50);
    }
    const seSpans = computeSESpans(seData);

    // ── 2. Extract R, G, B channels ───────────────────────────────────────────
    onProgress('2) Extracting colour channels…');
    await yield_();
    const R = new Float32Array(H * W);
    const G = new Float32Array(H * W);
    const B = new Float32Array(H * W);
    for (let i = 0; i < H * W; i++) {
        R[i] = rgba[i * 4];
        G[i] = rgba[i * 4 + 1];
        B[i] = rgba[i * 4 + 2];
    }

    // ── 3. Detect WSP cards (full resolution — mirrors Python build_card_mask) ─
    onProgress('3) Detecting WSP cards…');
    await yield_();

    onProgress('   Applying colour threshold…');
    await yield_();
    // mask = (R > B+50 & G > B+50) | (B > R+15 & B > G+15)  [create_card_mask.m]
    const rawMask = new Uint8Array(H * W);
    for (let i = 0; i < H * W; i++) {
        const r = R[i], g = G[i], b = B[i];
        rawMask[i] = ((r > b + 50 && g > b + 50) || (b > r + 15 && b > g + 15)) ? 1 : 0;
    }

    onProgress('   Filling holes…');
    await yield_();
    const filled = binaryFillHoles(rawMask, H, W);

    onProgress('   Removing small objects…');
    await yield_();
    const opened = removeSmallObjects(filled, H, W, MIN_CARD_AREA);

    onProgress('   Morphological closing (disk SE)…');
    await yield_();
    const cardMask = morphClose(opened, seSpans, H, W);

    onProgress('   Labelling card regions…');
    await yield_();
    const { labels, regions } = labelCC(cardMask, H, W);

    if (regions.length === 0)
        throw new Error('No card regions detected. Verify the image and DPI setting.');

    const cardRegions = orderCards(regions, W);
    onProgress(`   Found ${cardRegions.length} card(s):`);
    await yield_();
    for (let ci = 0; ci < cardRegions.length; ci++) {
        onProgress(`   • WSP-${ci + 1} detected`);
        await yield_();
    }

    // Greyscale image for card overlay visualisation
    const gray = new Uint8Array(H * W);
    for (let i = 0; i < H * W; i++)
        gray[i] = Math.round(0.299 * R[i] + 0.587 * G[i] + 0.114 * B[i]);

    // ── 4. Per-card analysis (full resolution) ────────────────────────────────
    const fullBinary = new Uint8Array(H * W);
    const cards      = [];

    for (let ci = 0; ci < cardRegions.length; ci++) {
        const cardRegion = cardRegions[ci];
        const cardLabel  = `WSP-${ci + 1}`;
        onProgress(`4.${ci + 1}) Analysing ${cardLabel} (${ci + 1}/${cardRegions.length})…`);
        await yield_();

        // Bounding box directly from full-resolution label map [mirrors Python card_region.bbox]
        const [rb0, cb0, rb1, cb1] = cardRegion.bbox;
        const cropH = rb1 - rb0, cropW = cb1 - cb0;

        // R crop + card mask from full-res labels [mirrors Python: mask_crop = (lbl == card_region.label)[rb0:rb1, cb0:cb1]]
        const R_crop    = new Float32Array(cropH * cropW);
        const mask_crop = new Uint8Array(cropH * cropW);
        for (let r = 0; r < cropH; r++)
            for (let c = 0; c < cropW; c++) {
                const fi = (rb0 + r) * W + (cb0 + c);
                R_crop[r * cropW + c]    = R[fi];
                mask_crop[r * cropW + c] = (labels[fi] === cardRegion.label) ? 1 : 0;
            }

        // binarize_card.m
        const bw = binarize(R_crop, mask_crop, cropH * cropW);

        // Accumulate full-image stain binary
        for (let r = 0; r < cropH; r++)
            for (let c = 0; c < cropW; c++)
                if (bw[r * cropW + c]) fullBinary[(rb0 + r) * W + (cb0 + c)] = 1;

        // Coverage
        let stainPx = 0, maskPx = 0;
        for (let i = 0; i < cropH * cropW; i++) { stainPx += bw[i]; maskPx += mask_crop[i]; }
        const coverage = maskPx > 0 ? 100.0 * stainPx / maskPx : 0;

        // detect_droplets.m + isolate_elements.m
        const initComps  = labelCC(bw, cropH, cropW).regions;
        const finalComps = isolateElements(initComps, R_crop, cropH, cropW);

        // droplet_statistics.m
        const statsRaw = dropletStatistics(finalComps, dpi, false);
        const statsSF  = dropletStatistics(finalComps, dpi, true);

        onProgress(`   ✓ ${cardLabel}: ${finalComps.length} droplets detected.`);
        await yield_();

        cards.push({
            label: cardLabel,
            bbox: [rb0, cb0, rb1, cb1],
            statsRaw, statsSF,
            coverage,
            nDroplets: finalComps.length,
            bw, mask_crop, cropH, cropW, rb0, cb0,
        });
    }

    return { H, W, cards, fullBinary, gray };
}
