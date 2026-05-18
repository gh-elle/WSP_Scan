function export_masks(image_path, output_folder)
%EXPORT_MASKS  Save per-card binary masks for direct comparison with Python.
%
%   Runs the same card-detection and binarisation pipeline as
%   compute_droplet_coverage.m, then writes four PNG files per card:
%
%     cardmask-WSP-N-<image>.png   — binary card region mask (0 / 255)
%     stainmask-WSP-N-<image>.png  — binary stain mask (R < 120 within card)
%     step1mask-WSP-N-<image>.png  — after detect_droplets  (labeled > 0;
%                                    should equal stainmask)
%     step2mask-WSP-N-<image>.png  — after isolate_elements (adaptive
%                                    sub-segmentation reconstruction)
%
%   step1mask and step2mask are diagnostic: comparing them with Python output
%   from --save-intermediate-masks reveals which pipeline stage causes any
%   MATLAB/Python numerical differences.
%
%   Note: binary_fill_holes (imfill) and bwareaopen are used only for
%   card-mask building (create_card_mask.m), not for stain analysis.
%   droplet_statistics applies no area filter, so step2mask corresponds to
%   the exact set of pixels used for droplet counting.
%
%   All images are cropped to the card's bounding box so they can be
%   compared pixel-by-pixel with Python output from:
%       python wsp_matlab.py <image> --save-masks --save-intermediate-masks
%                            [--results-folder <dir>]
%
%   Parameters
%   ----------
%   image_path    : full path to the scan image (e.g. '...20250718_A_UAS_56.tif')
%   output_folder : folder where PNGs are written (created if absent)
%
%   Example
%   -------
%   export_masks('C:\...\20250718_A_UAS_56.tif', 'C:\...\masks')

if nargin < 2
    error('Usage: export_masks(image_path, output_folder)');
end

if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

% ── Extract image name (no extension) ────────────────────────────────────────
[~, image_name, ~] = fileparts(image_path);

% ── Load image ────────────────────────────────────────────────────────────────
fprintf('Reading %s\n', image_path);
im = imread(image_path);
if size(im, 3) == 1
    im = repmat(im, 1, 1, 3);
end
R = im(:,:,1);
G = im(:,:,2);
B = im(:,:,3);

% ── Stage 1: card mask ────────────────────────────────────────────────────────
struct_elem = strel('disk', 50);
[card_masks, ~, card_labels] = create_card_mask(R, G, B, struct_elem, '', "");
close(gcf);    % close the figure opened by create_card_mask

n_cards = numel(card_masks);
fprintf('Found %d card(s)\n\n', n_cards);

% ── Stage 2: per-card binarisation and mask export ────────────────────────────
for idx = 1:n_cards
    lbl = card_labels{idx};
    slug = strrep(lbl, ' ', '-');   % e.g. "WSP-4"

    % Crop bounding box (same as Python: rows/cols where mask has any True pixel)
    row_mask = any(card_masks{idx}, 2);
    col_mask = any(card_masks{idx}, 1);

    card_crop  = card_masks{idx}(row_mask, col_mask);   % logical
    R_crop     = R(row_mask, col_mask);

    % binarize_card.m: stain = (R < 120) & card_mask
    stain_crop = (R_crop < 120) & card_crop;             % logical

    % Write binary PNGs (logical → uint8 0/255)
    cardmask_file  = fullfile(output_folder, sprintf('cardmask-%s-%s.png',  slug, image_name));
    stainmask_file = fullfile(output_folder, sprintf('stainmask-%s-%s.png', slug, image_name));

    imwrite(uint8(card_crop)  * 255, cardmask_file);
    imwrite(uint8(stain_crop) * 255, stainmask_file);

    % ── Intermediate masks: step1 (after detect_droplets) and step2 (after isolate_elements) ──

    % step1mask: detect_droplets on the stain binary — labeled > 0 should equal stainmask
    [~, ~, labeled_init] = detect_droplets(stain_crop);
    step1_crop = labeled_init > 0;
    step1mask_file = fullfile(output_folder, sprintf('step1mask-%s-%s.png', slug, image_name));
    imwrite(uint8(step1_crop) * 255, step1mask_file);

    % step2mask: isolate_elements — reconstruct binary from adaptive sub-segmentation.
    % For each initial connected component:
    %   • crop R_crop to the component bounding box
    %   • compute adaptive threshold = 190 - (190 - min(R_crop_bb)) / 2
    %   • run detect_droplets on (element < threshold)
    %   • if sub-components found: paint adaptive-threshold binary back to crop coords
    %   • if no sub-components: paint original stain pixels (component kept unchanged)
    [~, comp_init, ~] = detect_droplets(stain_crop);
    isol_crop = false(size(stain_crop));
    for k = 1:numel(comp_init)
        bb  = comp_init(k).BoundingBox;          % [xmin, ymin, width, height] (half-pixel offset)
        r0  = ceil(bb(2));                        % first row in card crop (1-indexed)
        c0  = ceil(bb(1));                        % first col in card crop (1-indexed)
        h   = round(bb(4));
        w   = round(bb(3));
        r1  = min(r0 + h - 1, size(isol_crop, 1));
        c1  = min(c0 + w - 1, size(isol_crop, 2));
        rh  = r1 - r0 + 1;   % actual height after boundary clamp
        cw  = c1 - c0 + 1;   % actual width  after boundary clamp

        elem    = double(R_crop(r0:r1, c0:c1));
        thresh  = 190 - (190 - min(elem(:))) / 2;
        sub_bw  = elem < thresh;
        [~, sub_comps, ~] = detect_droplets(sub_bw(1:rh, 1:cw));

        if numel(sub_comps) == 0
            % No sub-segmentation: keep original stain pixels
            isol_crop(r0:r1, c0:c1) = isol_crop(r0:r1, c0:c1) | stain_crop(r0:r1, c0:c1);
        else
            % Sub-segmented: use adaptive-threshold binary
            isol_crop(r0:r1, c0:c1) = isol_crop(r0:r1, c0:c1) | sub_bw(1:rh, 1:cw);
        end
    end
    step2mask_file = fullfile(output_folder, sprintf('step2mask-%s-%s.png', slug, image_name));
    imwrite(uint8(isol_crop) * 255, step2mask_file);

    fprintf('  %s: %d initial components -> %d initial detected  ->  %s\n', ...
        lbl, nnz(stain_crop), numel(comp_init), slug);
end

fprintf('\nMasks saved to: %s\n', output_folder);

end
