function [droplet_histogram, area_histogram, area_hist, ...
    diameters, mean_diameter, std_diameter, nmd, nmd10, nmd90, vmd, vmd10, vmd90, ...
    ch, csn, droplet_count, droplets_per_cm2] = ...
    droplet_statistics(components, diameter_limits, dpi, spread_factor)
%DROPLET_STATISTICS  Compute size and coverage statistics for detected droplets.
%
%   Parameters
%   ----------
%   components      : struct array returned by REGIONPROPS (Area field required)
%   diameter_limits : row vector of bin edges in microns, e.g. [50 100 150 200 300 400 500 600]
%   dpi             : scanner resolution in dots per inch (typically 600)
%   spread_factor   : logical — true applies the USDA-ARS spread-factor polynomial
%                     (same equation used by DepositScan) to convert stain diameters
%                     to real droplet diameters:
%                         droplet = 0.53549306 × stain − 0.000084839 × stain²
%
%   Returns (all diameters in microns)
%   -----------------------------------
%   droplet_histogram : n_bins × 1  droplet count per diameter bin
%   area_histogram    : n_bins × 1  total stain area per diameter bin (pixels)
%   area_hist         : 1 × max_area  droplet count indexed by area (pixels)
%   diameters         : 1 × max_area  equivalent diameter for each pixel-area value
%   mean_diameter, std_diameter : scalar mean and standard deviation (N-1)
%   nmd / nmd10 / nmd90         : numeric median diameter and percentiles
%   vmd / vmd10 / vmd90         : volume median diameter and percentiles
%   ch                          : coefficient of homogeneity (VMD / NMD)
%   csn                         : 1 × max_area  cumulative count (%)
%   droplet_count               : total number of detected droplets
%   droplets_per_cm2            : droplet density (card area = 2.6 × 7.6 cm)

pixel_length = 25.4 / (dpi / 1000); % pixel length in microns
card_area    = 2.6 * 7.6;           % card area in square cm

n_bins = numel(diameter_limits) + 1;

% Guard: no droplets detected on this card
% (~isfield handles the case where isolate_elements returns a bare struct())
if isempty(components) || ~isfield(components, 'Area')
    droplet_histogram = zeros(n_bins, 1);
    area_histogram    = zeros(n_bins, 1);
    area_hist         = [];
    diameters         = [];
    mean_diameter     = NaN;  std_diameter  = NaN;
    nmd  = NaN;  nmd10  = NaN;  nmd90  = NaN;
    vmd  = NaN;  vmd10  = NaN;  vmd90  = NaN;
    ch   = NaN;  csn    = [];
    droplet_count    = 0;
    droplets_per_cm2 = 0;
    return
end

% Build area-count histogram: area_hist(k) = number of droplets with area k pixels
areas     = [components.Area]';
area_hist = accumarray(areas, 1)';   % 1 × max(area) row vector

% Compute equivalent diameter for each pixel-area value (vectorised)
idx_vec    = 1:numel(area_hist);
stain_diam = pixel_length * sqrt(4 * idx_vec / pi);
if spread_factor
    % Apply USDA-ARS spread-factor correction (same equation used by DepositScan):
    %   droplet_diameter = 0.53549306 × stain_diameter − 0.000084839 × stain_diameter²
    % Reference: USDA Agricultural Research Service / DepositScan software.
    diameters = (0.53549306 * stain_diam) - (0.000084839 * stain_diam.^2);
else
    diameters = stain_diam;
end

% Expand to per-droplet diameter list and compute descriptive statistics
all_diameters = repelem(diameters, area_hist);
mean_diameter = mean(all_diameters);
std_diameter  = std(all_diameters);

% Bin every area index into its diameter range
% discretize uses [a,b) intervals so boundary values fall into the upper bin,
% matching the original loop logic (strict < for lower bin membership).
bin_idx = discretize(diameters, [-Inf, diameter_limits, Inf])';

droplet_histogram = accumarray(bin_idx, area_hist(:),              [n_bins, 1]);
area_histogram    = accumarray(bin_idx, area_hist(:) .* idx_vec',  [n_bins, 1]);

% Numeric Median Diameter (NMD) — index arithmetic preserved from original
cum_count = cumsum(area_hist);
total     = sum(area_hist);
csn       = 100 * cum_count / total;
nmd   = diameters(1 + sum(cum_count <= total / 2));
nmd10 = diameters(1 + sum(cum_count <= total / 10));
nmd90 = diameters(1 + sum(cum_count <= 9 * total / 10));

% Volume Median Diameter (VMD)
volumes   = area_hist .* ((4/3) * pi * (diameters / 2).^3);
cum_vol   = cumsum(volumes);
total_vol = sum(volumes);
vmd   = diameters(1 + sum(cum_vol <= total_vol / 2));
vmd10 = diameters(1 + sum(cum_vol <= total_vol / 10));
vmd90 = diameters(1 + sum(cum_vol <= 9 * total_vol / 10));

ch = vmd / nmd;

droplet_count    = sum(area_hist);
droplets_per_cm2 = droplet_count / card_area;

return
