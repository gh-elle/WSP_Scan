function [new_components, new_labeled] = isolate_elements(components, R, initial_threshold)
%ISOLATE_ELEMENTS  Adaptive sub-segmentation of merged droplet blobs.
%
%   For each detected component, crops the red channel to the component
%   bounding box and applies an adaptive threshold to separate overlapping
%   droplets.  If the tighter threshold produces no sub-regions the original
%   component is kept unchanged.
%
%   The adaptive threshold is:
%       threshold = initial_threshold - (initial_threshold - min(R_crop)) / 2
%
%   Parameters
%   ----------
%   components        : struct array from REGIONPROPS (BoundingBox and Area required)
%   R                 : red-channel crop of the current card (uint8 or double)
%   initial_threshold : starting R-channel threshold (recommended: 190)
%
%   Returns
%   -------
%   new_components : struct array of final (possibly sub-segmented) components
%   new_labeled    : cell array of labeled sub-images (one per original component)

fraction        = 2;
new_components  = cell(size(components));
new_labeled     = cell(size(components));

% For each detected element, try to sub-segment with a tighter threshold
component_count = 0;
for idx = 1:numel(components)
    element = imcrop(R, components(idx).BoundingBox);

    [~, new_components{idx}, new_labeled{idx}] = detect_droplets( ...
        element < initial_threshold - (initial_threshold - min(element(:))) / fraction);

    if numel(new_components{idx}) == 0
        new_components{idx} = components(idx);
        component_count     = component_count + 1;
    else
        component_count = component_count + numel(new_components{idx});
    end
end

nl    = cell(component_count, 1);
count = 1;
nc    = struct([]);   % empty struct array — avoids fieldless stub if no droplets exist
for idx = 1:numel(new_components)
    for idx2 = 1:numel(new_components{idx})
        if count == 1
            nc = new_components{idx}(idx2);
        else
            nc(count, 1) = new_components{idx}(idx2);
        end
        nl{count} = new_labeled{idx};
        count     = count + 1;
    end
end
new_components = nc;
new_labeled    = nl;

return
