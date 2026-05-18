function scan_cards(directory, results_dir, card_names)
%SCAN_CARDS  Process all WSP card images in a directory and compile statistics.
%
%   SCAN_CARDS(directory) processes every supported image file in 'directory'
%   and writes a summary CSV to '<directory>/results/'.
%
%   SCAN_CARDS(directory, results_dir) writes output to the specified folder.
%   Pass '' or [] to use the default location.
%
%   SCAN_CARDS(directory, results_dir, card_names) applies custom string labels
%   to each card in every image (e.g. {'A-H','B-H','C-H','A-L','B-L','C-L'}).
%
%   Supported image formats: .tif, .tiff, .png, .jpg, .jpeg, .bmp
%
%   Output
%   ------
%   <results_dir>/<directory_name>.csv — one row per card per image, columns:
%     filename, card,
%     NMD/NMD10/NMD90/VMD/VMD10/VMD90/CH/mean_diameter/std_diameter,
%     covered_area_pct/droplet_count/droplets_per_cm2,
%     and the same diameter columns with '_SF' suffix (spread-factor corrected).
%
%   <results_dir>/ also contains PNG plots for each image (see ANALYZE_CARD_IMAGE).
%
%   Example
%   -------
%     scan_cards('C:\scans\trial_01')
%     scan_cards('C:\scans\trial_01', 'C:\output\trial_01')
%     scan_cards('C:\scans\trial_01', '', {'A-H','B-H','C-H','A-L','B-L','C-L'})

% ── Defaults ──────────────────────────────────────────────────────────────────
if ~exist('results_dir', 'var') || isempty(results_dir)
    results_dir = fullfile(directory, 'results');
end
if ~exist('card_names', 'var')
    card_names = "";
end

% ── Output folder ─────────────────────────────────────────────────────────────
if ~exist(results_dir, 'dir')
    mkdir(results_dir);
end

% ── Column headers ─────────────────────────────────────────────────────────────
header = {'filename', 'card', ...
    'NMD', 'NMD10', 'NMD90', 'VMD', 'VMD10', 'VMD90', 'CH', ...
    'mean_diameter', 'std_diameter', ...
    'covered_area_pct', 'droplet_count', 'droplets_per_cm2', ...
    'NMD_SF', 'NMD10_SF', 'NMD90_SF', 'VMD_SF', 'VMD10_SF', 'VMD90_SF', ...
    'CH_SF', 'mean_diameter_SF', 'std_diameter_SF'};

% ── Supported image extensions ─────────────────────────────────────────────────
image_ext = {'.tif', '.tiff', '.png', '.jpg', '.jpeg', '.bmp'};

% ── Process images ────────────────────────────────────────────────────────────
file_list = dir(directory);
T = [];
for fi = 1:numel(file_list)
    fname = file_list(fi).name;
    [~, ~, ext] = fileparts(fname);
    if ~any(strcmpi(ext, image_ext))
        continue
    end

    [data_coverage, data_diameters, data_diameters_sf, card_labels] = ...
        analyze_card_image(directory, fname, results_dir, card_names);

    n_cards   = numel(card_labels);
    data_meta = cell(n_cards, 2);
    for ci = 1:n_cards
        data_meta{ci, 1} = fname;
        data_meta{ci, 2} = card_labels{ci};
    end

    row_table = cell2table([data_meta, data_diameters, data_coverage, data_diameters_sf]);
    T = [T; row_table]; %#ok<AGROW>
end

% ── Write CSV ─────────────────────────────────────────────────────────────────
if ~isempty(T)
    T.Properties.VariableNames = header;
    d = strtrim(directory);
    if ~isempty(d) && (d(end) == '\' || d(end) == '/')
        d(end) = [];
    end
    [~, folder_name] = fileparts(d);
    writetable(T, fullfile(results_dir, [folder_name, '.csv']));
end
