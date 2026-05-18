function export_subcomps_matlab(image_path, output_dir)
%EXPORT_SUBCOMPS_MATLAB  Write per-blob sub-component counts for all cards.
%
%   For each initial connected component on each card, runs the adaptive
%   sub-segmentation (same as isolate_elements.m) and writes a text file
%   analogous to the Python diagnostic export_subcomps_py.py.
%
%   Parameters
%   ----------
%   image_path : full path to the scan image
%   output_dir : folder where text files are written (one per card)
%
%   Example
%   -------
%   export_subcomps_matlab('C:\...\20250718_A_UAS_56.tif', 'C:\...\test')

if nargin < 2
    error('Usage: export_subcomps_matlab(image_path, output_dir)');
end

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

fprintf('Reading %s\n', image_path);
im = imread(image_path);
if size(im, 3) == 1
    im = repmat(im, 1, 1, 3);
end
R = im(:,:,1);
G = im(:,:,2);
B = im(:,:,3);

struct_elem = strel('disk', 50);
[card_masks, ~, card_labels] = create_card_mask(R, G, B, struct_elem, '', "");
close(gcf);

n_cards = numel(card_masks);
fprintf('Found %d card(s)\n\n', n_cards);

for idx = 1:n_cards
    lbl = card_labels{idx};
    slug = strrep(lbl, ' ', '-');

    row_mask = any(card_masks{idx}, 2);
    col_mask = any(card_masks{idx}, 1);

    bw = binarize_card(R, G, B, card_masks{idx}, 120);
    bw_crop = bw(row_mask, col_mask);
    R_crop  = R(row_mask, col_mask);

    [~, init_comps, ~] = detect_droplets(bw_crop);
    n_init = numel(init_comps);

    out_path = fullfile(output_dir, sprintf('matlab_subcomps_%s.txt', slug));
    fid = fopen(out_path, 'w');

    fprintf(fid, '# MATLAB sub-component diagnostic  card=%s\n', slug);
    fprintf(fid, '# n_initial=%d\n', n_init);
    fprintf(fid, '# idx  init_area  min_val    threshold  n_sub  sub_areas\n');

    total = 0;
    for k = 1:n_init
        % Use imcrop exactly as isolate_elements.m does
        elem   = double(imcrop(R_crop, init_comps(k).BoundingBox));
        min_v  = min(elem(:));
        thresh = 190 - (190 - min_v) / 2;
        sub_bw = elem < thresh;

        [~, sub_comps, ~] = detect_droplets(sub_bw);
        n_sub = numel(sub_comps);

        % Build sub_areas string
        if n_sub == 0
            areas_str = '[]';
        else
            areas = sort([sub_comps.Area], 'descend');
            areas_str = ['[' strjoin(arrayfun(@(x) num2str(x), areas, 'UniformOutput', false), ', ') ']'];
        end

        fprintf(fid, '  %4d  %10d  %9.1f  %9.4f  %5d  %s\n', ...
            k-1, init_comps(k).Area, min_v, thresh, n_sub, areas_str);

        if n_sub == 0
            total = total + 1;
        else
            total = total + n_sub;
        end
    end

    fprintf(fid, '# TOTAL=%d\n', total);
    fclose(fid);

    fprintf('  %s: %d initial -> %d final  saved to %s\n', slug, n_init, total, out_path);
end

end
