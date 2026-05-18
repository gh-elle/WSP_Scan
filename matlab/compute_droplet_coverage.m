function [nmd, nmd10, nmd90, vmd, vmd10, vmd90, ch, csn, diameters, mean_diameter, ...
    std_diameter, covered_area, droplet_count, droplets_per_cm2, ...
    droplet_histogram, card_labels, full_binary] = ...
    compute_droplet_coverage(filename, diameter_limits, plot_title, spread_factor, card_names)
%COMPUTE_DROPLET_COVERAGE  Core analysis loop: load image, detect droplets, compute stats.
%
%   Parameters
%   ----------
%   filename        : full path to the scan image (TIF/PNG/JPG …)
%   diameter_limits : bin edges in microns, e.g. [50 100 150 200 300 400 500 600]
%   plot_title      : string used as the card-segmentation figure title
%   spread_factor   : logical — true corrects stain diameters with the spread-factor polynomial
%   card_names      : string array of custom card labels, or "" for defaults ("WSP-1" …)
%
%   Scanner DPI is fixed at 600; edit the DROPLET_STATISTICS call to change it.
%
%   Returns  (all diameter values in microns)
%   -----------------------------------------
%   nmd … vmd90      : NMD/VMD family statistics (1 × n_cards)
%   ch               : coefficient of homogeneity VMD/NMD (1 × n_cards)
%   csn              : cumulative count (%) cell array (1 × n_cards)
%   diameters        : diameter vector cell array (1 × n_cards)
%   mean_diameter, std_diameter : per-card mean and std of droplet diameter (1 × n_cards)
%   covered_area     : percentage of card area covered by stain (1 × n_cards)
%   droplet_count    : total droplets per card (1 × n_cards)
%   droplets_per_cm2 : droplet density (1 × n_cards)
%   droplet_histogram: n_cards × n_bins droplet count histogram
%   card_labels      : 1 × n_cards cell of label strings
%   full_binary      : logical image — union of all card stain masks

fprintf('Reading %s\n', filename)
im = imread(filename);
if size(im, 3) == 1          % grayscale — replicate to 3 channels
    im = repmat(im, 1, 1, 3);
end
R  = im(:,:,1);
G  = im(:,:,2);
B  = im(:,:,3);

struct_elem = strel('disk', 50);
[card_masks, ~, card_labels] = create_card_mask(R, G, B, struct_elem, plot_title, card_names);

n_cards          = numel(card_masks);
mean_diameter    = zeros(1, n_cards);
std_diameter     = zeros(1, n_cards);
nmd              = zeros(1, n_cards);
nmd10            = zeros(1, n_cards);
nmd90            = zeros(1, n_cards);
vmd              = zeros(1, n_cards);
vmd10            = zeros(1, n_cards);
vmd90            = zeros(1, n_cards);
ch               = zeros(1, n_cards);
csn              = cell(1, n_cards);
diameters        = cell(1, n_cards);
covered_area     = zeros(1, n_cards);
droplet_count    = zeros(1, n_cards);
droplets_per_cm2 = zeros(1, n_cards);
droplet_histogram = zeros(n_cards, numel(diameter_limits) + 1);

fprintf('\nCard analysis:\t');
full_binary = false(size(B));
for idx = 1:n_cards
    fprintf('%s\t', card_labels{idx});

    % Pre-compute crop region once to avoid three repeated sum() calls
    row_mask = any(card_masks{idx}, 2);
    col_mask = any(card_masks{idx}, 1);

    bw = binarize_card(R, G, B, card_masks{idx}, 120);
    full_binary = full_binary | bw;

    [~, new_components, new_labeled] = detect_droplets(bw(row_mask, col_mask)); %#ok<ASGLU>

    [new_components, ~] = isolate_elements(new_components, R(row_mask, col_mask), 190);

    [droplet_histogram(idx,:), ~, ~, diameters{idx}, mean_diameter(idx), ...
        std_diameter(idx), nmd(idx), nmd10(idx), nmd90(idx), vmd(idx), ...
        vmd10(idx), vmd90(idx), ch(idx), csn{idx}, droplet_count(idx), ...
        droplets_per_cm2(idx)] = ...
        droplet_statistics(new_components, diameter_limits, 600, spread_factor);

    covered_area(idx) = 100 * sum(bw(:)) / sum(card_masks{idx}(:));
end
fprintf('\n');

return
