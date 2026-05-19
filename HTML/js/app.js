'use strict';
/**
 * app.js — UI controller for the WSP scan web app.
 * Handles: file upload, DPI selection, pipeline orchestration, chart rendering,
 * results table formatting, and CSV + PNG export.
 */

// ── DOM refs ──────────────────────────────────────────────────────────────────
const dropZone      = document.getElementById('drop-zone');
const fileInput     = document.getElementById('file-input');
const dpiSelect     = document.getElementById('dpi-select');
const customDpiRow  = document.getElementById('custom-dpi-row');
const customDpiInput= document.getElementById('custom-dpi-input');
const customDpiBtn  = document.getElementById('custom-dpi-btn');
const dpiWarning    = document.getElementById('dpi-warning');
const runBtn        = document.getElementById('run-btn');
const spinner       = document.getElementById('spinner');
const progressLog   = document.getElementById('progress-log');
const outputSection = document.getElementById('output-section');
const resultsTable  = document.getElementById('results-table');
const resultsContent= document.getElementById('results-content');
const exportSection = document.getElementById('export-section');
const exportBtn     = document.getElementById('export-btn');
const folderLabel   = document.getElementById('folder-label');
const pickFolderBtn = document.getElementById('pick-folder-btn');
const previewBox         = document.getElementById('preview-box');
const previewPlaceholder = document.getElementById('preview-placeholder');
const previewSpinner     = document.getElementById('preview-spinner');
const dropPreview        = document.getElementById('drop-preview');

// chart canvases
const cvOverlay   = document.getElementById('cv-overlay');
const cvStain     = document.getElementById('cv-stain');
const cvHist      = document.getElementById('cv-hist');
const cvCumul     = document.getElementById('cv-cumul');

// ── App state ─────────────────────────────────────────────────────────────────
let currentFile     = null;
let currentDpi      = 600;
let currentResult   = null;
let currentCharts   = {};
let folderHandle    = null;
let previewObjectURL= null;
let cachedImageData = null;  // decoded {rgba, H, W} — filled on file select for TIF

// ── FSA (File System Access API) availability ─────────────────────────────────
const fsaSupported = ('showDirectoryPicker' in window);
if (!fsaSupported) {
    pickFolderBtn.style.display = 'none';
    folderLabel.textContent = 'Files will download to your browser downloads folder';
}

// ── DPI selector ──────────────────────────────────────────────────────────────
function updateDpiWarning(dpi) {
    dpiWarning.style.display = dpi > 600 ? 'block' : 'none';
}

dpiSelect.addEventListener('change', () => {
    if (dpiSelect.value === 'custom') {
        customDpiRow.style.display = 'flex';
    } else {
        customDpiRow.style.display = 'none';
        currentDpi = parseInt(dpiSelect.value, 10);
        updateDpiWarning(currentDpi);
    }
});

customDpiBtn.addEventListener('click', () => {
    const v = parseInt(customDpiInput.value, 10);
    if (!v || v < 72 || v > 1200) {
        customDpiInput.style.borderColor = '#d32f2f';
        return;
    }
    customDpiInput.style.borderColor = '#43a047';
    currentDpi = v;
    updateDpiWarning(currentDpi);
});

// ── File input ────────────────────────────────────────────────────────────────
dropZone.addEventListener('click', () => fileInput.click());
fileInput.addEventListener('change', () => handleFile(fileInput.files[0]));

dropZone.addEventListener('dragover', e => { e.preventDefault(); dropZone.classList.add('drag-over'); });
dropZone.addEventListener('dragleave',  () => dropZone.classList.remove('drag-over'));
dropZone.addEventListener('drop', e => {
    e.preventDefault();
    dropZone.classList.remove('drag-over');
    handleFile(e.dataTransfer.files[0]);
});

async function handleFile(file) {
    if (!file) return;
    const ext = file.name.split('.').pop().toLowerCase();
    if (!['tif', 'tiff', 'png', 'jpg', 'jpeg'].includes(ext)) {
        alert('Unsupported format. Please choose a TIF, PNG, or JPEG file.');
        return;
    }
    currentFile     = file;
    cachedImageData = null;

    dropZone.classList.add('has-file');
    dropZone.querySelector('p').textContent     = file.name;
    dropZone.querySelector('.hint').textContent = 'Click or drop to change file';

    // Release previous object URL
    if (previewObjectURL) { URL.revokeObjectURL(previewObjectURL); previewObjectURL = null; }

    // Reset preview box state — make it visible now that a file is chosen
    previewBox.style.display         = 'flex';
    previewPlaceholder.style.display = 'none';
    previewSpinner.style.display     = 'none';
    dropPreview.style.display        = 'none';

    if (ext !== 'tif' && ext !== 'tiff') {
        // PNG / JPEG — browser renders natively, no decode needed
        previewObjectURL          = URL.createObjectURL(file);
        dropPreview.src           = previewObjectURL;
        dropPreview.style.display = 'block';
        runBtn.disabled = false;
    } else {
        // TIF — decode immediately and show thumbnail; enable button when done
        previewSpinner.style.display = 'flex';
        try {
            const { rgba, H, W } = await loadImage(file);
            cachedImageData = { rgba, H, W };
            previewSpinner.style.display = 'none';
            showTifPreview(rgba, H, W);
        } catch (e) {
            previewSpinner.style.display = 'none';
            previewPlaceholder.querySelector('p').textContent = '⚠ Could not decode preview';
            previewPlaceholder.style.display = 'flex';
        }
        runBtn.disabled = false;
    }
}

// ── Progress log helpers ──────────────────────────────────────────────────────
function logClear() {
    progressLog.innerHTML   = '';
    resultsTable.style.display  = 'none';
    resultsContent.style.display = 'none';
    exportSection.style.display = 'none';
    outputSection.style.display = '';   // show output section
}
function log(msg) {
    const line = document.createElement('div');
    line.textContent = msg;
    progressLog.appendChild(line);

    // Scroll window down so the bottom of the output section stays visible,
    // but never scrolls past it (stops when its bottom edge hits the viewport bottom).
    const rect = outputSection.getBoundingClientRect();
    if (rect.bottom > window.innerHeight) {
        window.scrollBy({ top: rect.bottom - window.innerHeight, behavior: 'smooth' });
    }
}

// ── Run analysis ──────────────────────────────────────────────────────────────
runBtn.addEventListener('click', runAnalysis);

async function runAnalysis() {
    if (!currentFile) return;
    if (dpiSelect.value === 'custom') {
        const v = parseInt(customDpiInput.value, 10);
        if (!v || v < 72 || v > 1200) { customDpiInput.style.borderColor = '#d32f2f'; return; }
        currentDpi = v;
    }

    runBtn.disabled       = true;
    spinner.style.display = 'block';
    logClear();
    log('0) Loading image…');

    try {
        // Reuse cached decode for TIF (already decoded during preview)
        const { rgba, H, W } = cachedImageData || await loadImage(currentFile);

        log(`   ✓ Image loaded: ${W} × ${H} px`);
        currentResult = await runPipeline(rgba, H, W, currentDpi, log);
        log('✓ Analysis complete.');
        displayResults(currentResult, currentDpi);
    } catch (err) {
        log(`✗ Error: ${err.message}`);
        console.error(err);
    } finally {
        runBtn.disabled       = false;
        spinner.style.display = 'none';
    }
}

/** Render a thumbnail from decoded RGBA into the preview box. */
function showTifPreview(rgba, H, W) {
    const MAX = 800;
    const scale = Math.min(1, MAX / Math.max(H, W));
    const tW = Math.round(W * scale), tH = Math.round(H * scale);
    const cv  = document.createElement('canvas');
    cv.width = tW; cv.height = tH;
    const ctx = cv.getContext('2d');
    const src = document.createElement('canvas');
    src.width = W; src.height = H;
    const srcCtx = src.getContext('2d');
    const id = srcCtx.createImageData(W, H);
    id.data.set(rgba.slice ? rgba.slice(0, W * H * 4) : rgba.subarray(0, W * H * 4));
    srcCtx.putImageData(id, 0, 0);
    ctx.drawImage(src, 0, 0, tW, tH);
    dropPreview.src           = cv.toDataURL('image/png');
    dropPreview.style.display = 'block';
}

// ── Image loader ──────────────────────────────────────────────────────────────
function loadImage(file) {
    return new Promise((resolve, reject) => {
        const ext = file.name.split('.').pop().toLowerCase();
        const reader = new FileReader();

        if (ext === 'tif' || ext === 'tiff') {
            reader.onload = e => {
                if (typeof UTIF === 'undefined') {
                    reject(new Error('UTIF library not loaded. An internet connection is required for TIF support.'));
                    return;
                }
                const buf  = e.target.result;
                const ifds = UTIF.decode(buf);
                UTIF.decodeImage(buf, ifds[0]);
                const ifd  = ifds[0];
                const W    = ifd.width, H = ifd.height;
                const rgba = UTIF.toRGBA8(ifd);
                resolve({ rgba, H, W });
            };
            reader.onerror = () => reject(new Error('Failed to read file.'));
            reader.readAsArrayBuffer(file);
        } else {
            reader.onload = e => {
                const img = new Image();
                img.onload = () => {
                    const cv  = document.createElement('canvas');
                    cv.width  = img.width; cv.height = img.height;
                    const ctx = cv.getContext('2d');
                    ctx.drawImage(img, 0, 0);
                    const data = ctx.getImageData(0, 0, img.width, img.height);
                    resolve({ rgba: data.data, H: img.height, W: img.width });
                };
                img.onerror = () => reject(new Error('Failed to decode image.'));
                img.src = e.target.result;
            };
            reader.onerror = () => reject(new Error('Failed to read file.'));
            reader.readAsDataURL(file);
        }
    });
}

// ── Display results ───────────────────────────────────────────────────────────
function displayResults(result, dpi) {
    resultsTable.textContent     = formatTable(result.cards);
    resultsTable.style.display   = 'block';
    resultsContent.style.display = '';
    exportSection.style.display  = '';

    destroyCharts();
    renderOverlay(result);
    renderStainMask(result);
    renderHistogram(result.cards);
    renderCumulative(result.cards);

    // Scroll so the bottom of the terminal block (including the results table)
    // aligns with the bottom of the viewport — wait one frame for layout to settle.
    requestAnimationFrame(() => {
        const rect = outputSection.getBoundingClientRect();
        if (rect.bottom > window.innerHeight) {
            window.scrollBy({ top: rect.bottom - window.innerHeight, behavior: 'smooth' });
        }
    });
}

// ── Results table (matches Python stdout exactly) ─────────────────────────────
function formatTable(cards) {
    const LABEL_W = 22, VAL_W = 9;
    const n = cards.length;
    const hdr = 'Card analysis:'.padEnd(LABEL_W) +
                cards.map(c => c.label.padStart(VAL_W)).join('');
    const sep = '-'.repeat(LABEL_W + VAL_W * n);

    const row = (label, fn) =>
        label.padEnd(LABEL_W) + cards.map(c => fn(c).padStart(VAL_W)).join('');

    const fv  = v  => (v == null || isNaN(v)) ? '   N/A' : v.toFixed(2);
    const fi  = v  => String(Math.round(v));

    return [
        hdr, sep,
        'Without Spread Factor',
        row('NMD',          c => fv(c.statsRaw?.nmd)),
        row('NMD10',        c => fv(c.statsRaw?.nmd10)),
        row('NMD90',        c => fv(c.statsRaw?.nmd90)),
        row('VMD',          c => fv(c.statsRaw?.vmd)),
        row('VMD10',        c => fv(c.statsRaw?.vmd10)),
        row('VMD90',        c => fv(c.statsRaw?.vmd90)),
        row('CH',           c => fv(c.statsRaw?.ch)),
        row('Mean diameter',c => fv(c.statsRaw?.mean)),
        row('Diameter STD', c => fv(c.statsRaw?.std)),
        sep,
        'With Spread Factor',
        row('NMD',          c => fv(c.statsSF?.nmd)),
        row('NMD10',        c => fv(c.statsSF?.nmd10)),
        row('NMD90',        c => fv(c.statsSF?.nmd90)),
        row('VMD',          c => fv(c.statsSF?.vmd)),
        row('VMD10',        c => fv(c.statsSF?.vmd10)),
        row('VMD90',        c => fv(c.statsSF?.vmd90)),
        row('CH',           c => fv(c.statsSF?.ch)),
        row('Mean diameter',c => fv(c.statsSF?.mean)),
        row('Diameter STD', c => fv(c.statsSF?.std)),
        sep,
        row('Covered area',   c => c.coverage.toFixed(2) + '%'),
        row('Droplets per card',   c => fi(c.nDroplets)),
        row('Droplets per cm2',    c => fv(c.nDroplets / CARD_AREA_CM2)),
        sep,
    ].join('\n');
}

// ── Chart helpers ─────────────────────────────────────────────────────────────
function destroyCharts() {
    for (const ch of Object.values(currentCharts)) ch.destroy();
    currentCharts = {};
}

/** Card overlay: greyscale image (display-scaled) with black bounding boxes and WSP labels. */
function renderOverlay(result) {
    const { H, W, cards, gray } = result;

    // Scale to same display size as the stain mask (~2 MP max)
    const scale = Math.min(1, Math.sqrt(2e6 / (H * W)));
    const cH = Math.round(H * scale), cW = Math.round(W * scale);
    cvOverlay.width  = cW;
    cvOverlay.height = cH;
    const ctx = cvOverlay.getContext('2d');

    // Draw greyscale background (nearest-neighbour downsample)
    const imgData = ctx.createImageData(cW, cH);
    const d = imgData.data;
    for (let r = 0; r < cH; r++) {
        const sr = Math.floor(r / scale);
        for (let c = 0; c < cW; c++) {
            const sc = Math.floor(c / scale);
            const g  = gray[sr * W + sc];
            const oi = (r * cW + c) * 4;
            d[oi] = g; d[oi+1] = g; d[oi+2] = g; d[oi+3] = 255;
        }
    }
    ctx.putImageData(imgData, 0, 0);

    // Draw bounding box + label for each card
    const fontSize  = Math.max(10, Math.round(cW / 55));
    const lineW     = Math.max(1,  Math.round(cW / 350));
    ctx.font         = `bold ${fontSize}px sans-serif`;
    ctx.textAlign    = 'center';

    cards.forEach((card) => {
        const [r0f, c0f, r1f, c1f] = card.bbox;
        // Scale full-res bbox to display canvas
        const r0 = Math.round(r0f * scale), c0 = Math.round(c0f * scale);
        const r1 = Math.round(r1f * scale), c1 = Math.round(c1f * scale);
        // Bounding box
        ctx.strokeStyle = '#000000';
        ctx.lineWidth   = lineW;
        ctx.strokeRect(c0, r0, c1 - c0, r1 - r0);
        // Label above the top edge with white halo for legibility
        const cx = (c0 + c1) / 2;
        const ty = r0 - lineW - 2;
        ctx.textBaseline = 'bottom';
        ctx.lineWidth   = fontSize * 0.35;
        ctx.strokeStyle = '#ffffff';
        ctx.strokeText(card.label, cx, ty);
        ctx.fillStyle   = '#000000';
        ctx.fillText(card.label, cx, ty);
    });
}

/** Stain mask: black stains on white background with thin red card segmentation borders. */
function renderStainMask(result) {
    const { H, W, fullBinary, cards } = result;
    const scale = Math.min(1, Math.sqrt(2e6 / (H * W)));
    const cH = Math.round(H * scale), cW = Math.round(W * scale);
    cvStain.width  = cW;
    cvStain.height = cH;
    const ctx = cvStain.getContext('2d');
    const img = ctx.createImageData(cW, cH);
    const d   = img.data;

    // Helper: is canvas pixel (r,c) inside card's mask_crop?
    const insideMask = (card, r, c) => {
        const sr = Math.floor(r / scale), sc = Math.floor(c / scale);
        const lr = sr - card.rb0, lc = sc - card.cb0;
        if (lr < 0 || lr >= card.cropH || lc < 0 || lc >= card.cropW) return false;
        return card.mask_crop[lr * card.cropW + lc] === 1;
    };

    const DN = [[-1, 0], [1, 0], [0, -1], [0, 1]];

    for (let r = 0; r < cH; r++) {
        const sr = Math.floor(r / scale);
        for (let c = 0; c < cW; c++) {
            const sc = Math.floor(c / scale);
            const oi = (r * cW + c) * 4;

            // Border: canvas pixel is inside a mask AND a canvas neighbour is outside
            let isBorder = false;
            for (const card of cards) {
                if (!insideMask(card, r, c)) continue;
                for (const [dr, dc] of DN) {
                    const nr = r + dr, nc = c + dc;
                    if (nr < 0 || nr >= cH || nc < 0 || nc >= cW || !insideMask(card, nr, nc)) {
                        isBorder = true;
                        break;
                    }
                }
                break;
            }

            if (isBorder) {
                d[oi] = 255; d[oi+1] = 0; d[oi+2] = 0; d[oi+3] = 255;   // red
            } else {
                const g = fullBinary[sr * W + sc] ? 0 : 255;  // stain=black, bg=white
                d[oi] = g; d[oi+1] = g; d[oi+2] = g; d[oi+3] = 255;
            }
        }
    }
    ctx.putImageData(img, 0, 0);
}

/** Histogram: grouped bar chart per card (NMD bins). */
function renderHistogram(cards) {
    const labels = ['<50', '50–100', '100–150', '150–200', '200–300', '300–400', '400–500', '500–600', '>600'];
    const datasets = cards.map((card, i) => ({
        label:           card.label,
        data:            Array.from(card.statsRaw?.dropletHistogram ?? new Int32Array(9)),
        backgroundColor: hexWithAlpha(COLORS_HEX[i % COLORS_HEX.length], 0.7),
        borderColor:     COLORS_HEX[i % COLORS_HEX.length],
        borderWidth: 1,
    }));
    currentCharts.hist = new Chart(cvHist, {
        type: 'bar',
        data: { labels, datasets },
        options: {
            responsive: true,
            maintainAspectRatio: true,
            aspectRatio: 8 / 6,
            plugins: { legend: { position: 'top', labels: { font: { size: 11 } } } },
            scales: {
                x: { title: { display: true, text: 'Stain diameter (μm)' } },
                y: { title: { display: true, text: 'Droplet count' }, beginAtZero: true },
            },
        },
    });
}

/** Cumulative distribution: scatter+line for each card (spread-factor). */
function renderCumulative(cards) {
    const datasets = [];
    for (let i = 0; i < cards.length; i++) {
        const card = cards[i];
        const sf   = card.statsSF;
        if (!sf) continue;
        const points = [];
        for (let k = 0; k < sf.diam.length; k++)
            if (sf.hist[k] > 0) points.push({ x: sf.diam[k], y: sf.csn[k] });
        datasets.push({
            label:       card.label,
            data:        points,
            borderColor: COLORS_HEX[i % COLORS_HEX.length],
            backgroundColor: 'transparent',
            showLine: true,
            tension:  0.1,
            pointRadius: 0,
            borderWidth: 1.8,
        });
    }
    currentCharts.cumul = new Chart(cvCumul, {
        type: 'scatter',
        data: { datasets },
        options: {
            responsive: true,
            maintainAspectRatio: true,
            aspectRatio: 8 / 6,
            parsing: false,
            plugins: { legend: { position: 'top', labels: { font: { size: 11 } } } },
            scales: {
                x: { type: 'linear', title: { display: true, text: 'Diameter with SF (μm)' } },
                y: { title: { display: true, text: 'Cumulative count (%)' }, min: 0, max: 100 },
            },
        },
    });
}

// ── Export ────────────────────────────────────────────────────────────────────
pickFolderBtn.addEventListener('click', async () => {
    if (!fsaSupported) return;
    try {
        folderHandle = await window.showDirectoryPicker({ mode: 'readwrite' });
        folderLabel.textContent = '✓ ' + folderHandle.name;
    } catch { /* user cancelled */ }
});

exportBtn.addEventListener('click', exportResults);

async function exportResults() {
    if (!currentResult) return;

    // Try to open folder picker if not already set.
    // FSA is only available in Chrome/Edge; Firefox falls back to individual downloads.
    if (!folderHandle && fsaSupported) {
        folderLabel.textContent = 'Navigate to the html/results folder…';
        try {
            folderHandle = await window.showDirectoryPicker({ mode: 'readwrite' });
            folderLabel.textContent = '✓ ' + folderHandle.name;
        } catch {
            // Picker blocked or cancelled — fall through to individual downloads
            folderLabel.textContent = 'Files will download to your browser downloads folder';
        }
    }
    exportBtn.disabled = true;
    try {
        const csvData  = buildCSV(currentResult.cards, currentDpi);
        const txtData  = formatTable(currentResult.cards);
        const doCSV    = document.getElementById('chk-csv').checked;
        const doPlots  = document.getElementById('chk-plots').checked;

        if (doCSV) {
            await saveFile('wsp_results.csv', new Blob([csvData], { type: 'text/csv' }));
            await saveFile('wsp_results.txt', new Blob([txtData], { type: 'text/plain' }));
        }
        if (doPlots) {
            await saveCanvasAsPng(cvOverlay, 'plot_card_overlay.png');
            await saveCanvasAsPng(cvStain,   'plot_stain_mask.png');
            if (currentCharts.hist)  await saveChartAsPng(currentCharts.hist,  'plot_histogram.png');
            if (currentCharts.cumul) await saveChartAsPng(currentCharts.cumul, 'plot_cumulative.png');
        }
    } finally {
        exportBtn.disabled = false;
    }
}

function buildCSV(cards, dpi) {
    const rows = [
        ['Card', 'DPI',
         'NMD_raw','NMD10_raw','NMD90_raw','VMD_raw','VMD10_raw','VMD90_raw','CH_raw','Mean_raw','Std_raw',
         'NMD_sf', 'NMD10_sf', 'NMD90_sf', 'VMD_sf', 'VMD10_sf', 'VMD90_sf', 'CH_sf', 'Mean_sf', 'Std_sf',
         'Coverage_pct','Droplets','Droplets_per_cm2'].join(','),
    ];
    for (const c of cards) {
        const r = c.statsRaw, s = c.statsSF;
        const f = v => v == null || isNaN(v) ? '' : v.toFixed(4);
        rows.push([
            c.label, dpi,
            f(r?.nmd),f(r?.nmd10),f(r?.nmd90),f(r?.vmd),f(r?.vmd10),f(r?.vmd90),f(r?.ch),f(r?.mean),f(r?.std),
            f(s?.nmd),f(s?.nmd10),f(s?.nmd90),f(s?.vmd),f(s?.vmd10),f(s?.vmd90),f(s?.ch),f(s?.mean),f(s?.std),
            c.coverage.toFixed(4), c.nDroplets, (c.nDroplets / CARD_AREA_CM2).toFixed(4),
        ].join(','));
    }
    return rows.join('\n');
}

async function saveFile(name, blob) {
    if (folderHandle) {
        try {
            const fh = await folderHandle.getFileHandle(name, { create: true });
            const wr = await fh.createWritable();
            await wr.write(blob);
            await wr.close();
            return;
        } catch { /* fall through to download */ }
    }
    const url = URL.createObjectURL(blob);
    const a   = document.createElement('a');
    a.href    = url; a.download = name; a.click();
    setTimeout(() => URL.revokeObjectURL(url), 10000);
}

async function saveCanvasAsPng(canvas, name) {
    const blob = await new Promise(res => canvas.toBlob(res, 'image/png'));
    if (blob) await saveFile(name, blob);
}

async function saveChartAsPng(chart, name) {
    await saveCanvasAsPng(chart.canvas, name);
}

// ── Colour utilities ──────────────────────────────────────────────────────────
function hexToRGB(hex) {
    const v = parseInt(hex.slice(1), 16);
    return [(v >> 16) & 0xff, (v >> 8) & 0xff, v & 0xff];
}

function hexWithAlpha(hex, a) {
    const [r, g, b] = hexToRGB(hex);
    return `rgba(${r},${g},${b},${a})`;
}
