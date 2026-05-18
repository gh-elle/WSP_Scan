function dump_blob_matlab(image_path, output_dir)
%DUMP_BLOB_MATLAB  Dump imcrop and direct-index crops for the discrepant blob.
%
%  Finds the WSP-2 initial component with Area==9 and min R==36,
%  then prints both the imcrop result and the direct R_crop(r0:r1,c0:c1)
%  result so they can be compared against Python's bbox slice.
%
%  Usage:
%    dump_blob_matlab('C:\...\20250718_A_UAS_56.tif', 'C:\...\test')

if nargin < 2, error('Usage: dump_blob_matlab(image_path, output_dir)'); end

im = imread(image_path);
R  = im(:,:,1);
G  = im(:,:,2);
B  = im(:,:,3);

struct_elem = strel('disk', 50);
[card_masks, ~, card_labels] = create_card_mask(R, G, B, struct_elem, '', "");
close(gcf);

% Find WSP-2
wsp2 = 0;
for i = 1:numel(card_labels)
    if strcmp(card_labels{i}, 'WSP-2'), wsp2 = i; break; end
end
if wsp2 == 0, error('WSP-2 not found'); end

row_mask = any(card_masks{wsp2}, 2);
col_mask = any(card_masks{wsp2}, 1);

bw     = binarize_card(R, G, B, card_masks{wsp2}, 120);
bw_crop = bw(row_mask, col_mask);
R_crop  = R(row_mask, col_mask);

[~, init_comps, ~] = detect_droplets(bw_crop);

for k = 1:numel(init_comps)
    if init_comps(k).Area ~= 9, continue; end
    bb   = init_comps(k).BoundingBox;
    elem_imcrop = double(imcrop(R_crop, bb));
    if min(elem_imcrop(:)) ~= 36.0, continue; end

    % Direct indexing (old method)
    r0 = ceil(bb(2)); c0 = ceil(bb(1));
    h  = round(bb(4)); w  = round(bb(3));
    r1 = min(r0+h-1, size(R_crop,1));
    c1 = min(c0+w-1, size(R_crop,2));
    elem_direct = double(R_crop(r0:r1, c0:c1));

    thresh = 190 - (190 - 36.0) / 2;   % = 113

    fprintf('\n=== Discrepant blob (MATLAB idx %d, Area=9, min=36) ===\n', k-1);
    fprintf('BoundingBox = [%.1f %.1f %.1f %.1f]  (x y w h = col row ncols nrows)\n', ...
            bb(1), bb(2), bb(3), bb(4));
    fprintf('\nimcrop shape: %d x %d\n', size(elem_imcrop,1), size(elem_imcrop,2));
    disp(elem_imcrop);
    fprintf('direct index shape: %d x %d\n', size(elem_direct,1), size(elem_direct,2));
    disp(elem_direct);

    sub_bw_imcrop  = elem_imcrop  < thresh;
    sub_bw_direct  = elem_direct  < thresh;
    fprintf('sub_bw_imcrop (< %g):\n', thresh); disp(sub_bw_imcrop);
    fprintf('sub_bw_direct (< %g):\n', thresh); disp(sub_bw_direct);

    [~, c_ic, ~] = detect_droplets(sub_bw_imcrop);
    [~, c_di, ~] = detect_droplets(sub_bw_direct);
    fprintf('n_sub via imcrop       = %d\n', numel(c_ic));
    fprintf('n_sub via direct index = %d\n', numel(c_di));

    csvwrite(fullfile(output_dir, 'blob_crop_matlab_imcrop.csv'),  elem_imcrop);
    csvwrite(fullfile(output_dir, 'blob_crop_matlab_direct.csv'),  elem_direct);
    fprintf('\nSaved CSVs to %s\n', output_dir);
    break
end

end
