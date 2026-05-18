function [card_masks, components, card_labels] = ...
    create_card_mask(R, G, B, struct_elem, plot_title, card_names)
%CREATE_CARD_MASK  Segment WSP cards from the scan image using RGB thresholding.
%
%   The two detection conditions are:
%     R > B+50  AND  G > B+50   — yellow unsprayed card background
%     B > R+15  AND  B > G+15   — blue/purple water-stained card regions
%   Both tests use MATLAB's saturating uint8 integer arithmetic.
%   Cards are separated into a left column and a right column by position,
%   then sorted top-to-bottom within each column.
%
%   Parameters
%   ----------
%   R, G, B      : uint8 channel matrices (size H × W) from IMREAD
%   struct_elem  : morphological structuring element for binary closing (STREL)
%   plot_title   : string shown as the segmentation figure title
%   card_names   : string array of custom labels, or "" for defaults ("WSP-1" …)
%
%   Returns
%   -------
%   card_masks  : 1 × n_cards cell of logical masks (one per card, size H × W)
%   components  : 1 × n_cards struct array from REGIONPROPS
%   card_labels : 1 × n_cards cell of label strings

card_mask = (R > B + 50 & G > B + 50) | (B > R + 15 & B > G + 15);
card_mask = imfill(card_mask, 'holes');

% Remove small objects from binary image
card_mask = bwareaopen(card_mask, 1000);
card_mask = imclose(card_mask, struct_elem);

cc         = bwconncomp(card_mask, 8);
components = regionprops(cc, 'all');
labeled    = labelmatrix(cc);

n          = numel(components);
card_masks = cell(1, n);
x_pos      = zeros(1, n);
y_pos      = zeros(1, n);
for idx = 1:n
    card_masks{idx} = labeled == idx;
    x_pos(idx) = components(idx).BoundingBox(1);
    y_pos(idx) = components(idx).BoundingBox(2);
end

% Split components into left and right columns by position, then sort each
% column top-to-bottom. Using find() avoids assuming regionprops order.
left_range = fix(size(card_mask, 2) / 3);
left_idx   = find(x_pos <  left_range);
right_idx  = find(x_pos >= left_range);
[~, sort_left]  = sort(y_pos(left_idx));
[~, sort_right] = sort(y_pos(right_idx));
idx = [left_idx(sort_left), right_idx(sort_right)];

% Pre-allocate a cell array for card labels
card_labels = cell(1, n);
if card_names == ""
    for i = 1:n
        card_labels{i} = sprintf('WSP-%d', i);
    end
else
    for i = 1:n
        card_labels{i} = sprintf('WSP %s', card_names(i));
    end
end

card_masks = card_masks(idx);
components = components(idx);

imshow(card_mask)
for i = 1:numel(components)
    coord = int16(components(i).Centroid);
    text(coord(1), coord(2), sprintf('%s', card_labels{i}), ...
        'Color', 'red', 'FontSize', 14, 'HorizontalAlignment', 'center')
end
title(plot_title);

return
