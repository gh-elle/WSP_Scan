function [cc, components, labeled] = detect_droplets(bw)
%DETECT_DROPLETS  Label connected components in a binary stain image.
%
%   Uses 8-connectivity (MATLAB BWCONNCOMP default for 2-D images).
%   No area filter is applied; all components are returned regardless of size.
%
%   Parameters
%   ----------
%   bw : logical binary image (true = stain pixel)
%
%   Returns
%   -------
%   cc         : connected-components struct from BWCONNCOMP
%   components : struct array from REGIONPROPS (all properties)
%   labeled    : label matrix from LABELMATRIX

cc         = bwconncomp(bw, 8);
labeled    = labelmatrix(cc);
components = regionprops(cc, 'all');

return
