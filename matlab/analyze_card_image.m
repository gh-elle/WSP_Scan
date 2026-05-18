function [data_coverage, data_diameters, data_diameters_sf, card_labels] = ...
    analyze_card_image(folder, filename, results_folder, card_names)
%ANALYZE_CARD_IMAGE  Run the full WSP pipeline on a single scan image.
%
%   Calls COMPUTE_DROPLET_COVERAGE twice: once for raw stain diameters and
%   once with the spread-factor correction.  Saves four PNG files to
%   results_folder: card-segmentation overlay, binary stain mask, droplet
%   count histogram, and cumulative diameter distribution.
%
%   Parameters
%   ----------
%   folder         : directory containing the image file
%   filename       : image filename (e.g. 'scan_001.tif')
%   results_folder : output folder for PNG plots
%   card_names     : string array of custom card labels, or "" for defaults
%
%   Returns
%   -------
%   data_coverage     : n_cards × 3 cell  {covered_area_pct, droplet_count, droplets_per_cm2}
%   data_diameters    : n_cards × 9 cell  {NMD,NMD10,NMD90,VMD,VMD10,VMD90,CH,mean,std}
%   data_diameters_sf : same as data_diameters, with spread-factor correction
%   card_labels       : 1 × n_cards cell of label strings

name_parts = split(filename, '.');
name_parts = split(name_parts{1}, '_');
name_parts = join(name_parts, '-');
image_name = name_parts{1};

filepath = fullfile(folder, filename);

diameter_limits = [50, 100, 150, 200, 300, 400, 500, 600];
categories = cell(numel(diameter_limits) + 1, 1);
categories{1} = sprintf('<%d', diameter_limits(1));
for idx = 1:numel(diameter_limits) - 1
    categories{idx+1} = sprintf('%d-%d', diameter_limits(idx), diameter_limits(idx+1));
end
categories{end} = sprintf('>%d', diameter_limits(end));
diam_cat = categorical(categories);
diam_cat = reordercats(diam_cat, categories);

colors = {[1,0,0],[0,0,1],[0,1,0],[0,1,1],[1,0,1],[1,0.5,0],[0.5,0.5,1], ...
          [1,0.5,0.5],[1,1,0],[0.5,0.5,1],[0.5,1,0.5],[0.3,1,0.3],[0.5,0.5,0.5]};

[nmd, nmd10, nmd90, vmd, vmd10, vmd90, ch, csn, diameters, mean_diameter, std_diameter, ...
    covered_area, droplet_count, droplets_per_cm2, ...
    droplet_histogram, card_labels, full_binary] = ...
    compute_droplet_coverage(filepath, diameter_limits, image_name, false, card_names);

[nmd_sf, nmd10_sf, nmd90_sf, vmd_sf, vmd10_sf, vmd90_sf, ch_sf, csn_sf, diameters_sf, ...
    mean_diameter_sf, std_diameter_sf, covered_area, droplet_count, droplets_per_cm2, ...
    droplet_histogram, card_labels, full_binary] = ...
    compute_droplet_coverage(filepath, diameter_limits, image_name, true, card_names);

saveas(gcf, fullfile(results_folder, ['cards-', image_name, '.png']));
imwrite(full_binary, fullfile(results_folder, ['mask-cards-', image_name, '.png']));

figure
bar_handle = bar(diam_cat, droplet_histogram');
for bi = 1:size(bar_handle, 2)
    bar_handle(bi).FaceColor = colors{bi};
end
legend_labels = card_labels;
title('Droplets Number');
xlabel('Diameter (\mum)')
ylabel('Droplets')
legend(legend_labels, 'Location', 'northeastoutside')
exportgraphics(gcf, fullfile(results_folder, ['hist-', image_name, '.png']), 'Resolution', 300)

max_diam = 1200;
x_limit  = 0;
for ci = 1:size(card_labels, 2)
    x_limit = max(x_limit, sum(diameters{ci} < max_diam));
end
figure
hold on
for ci = 1:size(card_labels, 2)
    if x_limit > size(diameters{ci}, 2)
        plot([0, diameters{ci}], [0, csn{ci}], 'color', colors{ci}, 'LineWidth', 1)
    else
        plot([0, diameters{ci}(1:x_limit)], [0, csn{ci}(1:x_limit)], ...
            'color', colors{ci}, 'LineWidth', 1)
    end
end
title('Cumulative Sum');
xlabel('Diameter (\mum)')
ylabel('Droplets (%)')
xticks(0:100:max_diam)
legend(legend_labels, 'Location', 'southeast')
exportgraphics(gcf, fullfile(results_folder, ['cumulative-', image_name, '.png']), 'Resolution', 300)

fprintf('Without Spread Factor');
fprintf('\nNMD');            fprintf('\t%4.2f', nmd);
fprintf('\nNMD10');          fprintf('\t%4.2f', nmd10);
fprintf('\nNMD90');          fprintf('\t%4.2f', nmd90);
fprintf('\nVMD');            fprintf('\t%4.2f', vmd);
fprintf('\nVMD10');          fprintf('\t%4.2f', vmd10);
fprintf('\nVMD90');          fprintf('\t%4.2f', vmd90);
fprintf('\nCH');             fprintf('\t%4.2f', ch);
fprintf('\nMean diameter');  fprintf('\t%4.2f', mean_diameter);
fprintf('\nDiameter STD');   fprintf('\t%4.2f', std_diameter);

fprintf('\nWith Spread Factor');
fprintf('\nNMD');            fprintf('\t%4.2f', nmd_sf);
fprintf('\nNMD10');          fprintf('\t%4.2f', nmd10_sf);
fprintf('\nNMD90');          fprintf('\t%4.2f', nmd90_sf);
fprintf('\nVMD');            fprintf('\t%4.2f', vmd_sf);
fprintf('\nVMD10');          fprintf('\t%4.2f', vmd10_sf);
fprintf('\nVMD90');          fprintf('\t%4.2f', vmd90_sf);
fprintf('\nCH');             fprintf('\t%4.2f', ch_sf);
fprintf('\nMean diameter');  fprintf('\t%4.2f', mean_diameter_sf);
fprintf('\nDiameter STD');   fprintf('\t%4.2f', std_diameter_sf);

fprintf('\nCovered area');       fprintf('\t%4.2f%%', covered_area);
fprintf('\nDroplets per card');  fprintf('\t%4.2d',   droplet_count);
fprintf('\nDroplets per cm2');   fprintf('\t%4.2f',   droplets_per_cm2);
fprintf('\n');

n_cards          = numel(card_labels);
data_coverage    = cell(n_cards, 3);
data_diameters   = cell(n_cards, 9);
data_diameters_sf = cell(n_cards, 9);
for idx = 1:n_cards
    data_coverage{idx, 1} = covered_area(idx);
    data_coverage{idx, 2} = droplet_count(idx);
    data_coverage{idx, 3} = droplets_per_cm2(idx);
    data_diameters{idx, 1} = nmd(idx);
    data_diameters{idx, 2} = nmd10(idx);
    data_diameters{idx, 3} = nmd90(idx);
    data_diameters{idx, 4} = vmd(idx);
    data_diameters{idx, 5} = vmd10(idx);
    data_diameters{idx, 6} = vmd90(idx);
    data_diameters{idx, 7} = ch(idx);
    data_diameters{idx, 8} = mean_diameter(idx);
    data_diameters{idx, 9} = std_diameter(idx);
    data_diameters_sf{idx, 1} = nmd_sf(idx);
    data_diameters_sf{idx, 2} = nmd10_sf(idx);
    data_diameters_sf{idx, 3} = nmd90_sf(idx);
    data_diameters_sf{idx, 4} = vmd_sf(idx);
    data_diameters_sf{idx, 5} = vmd10_sf(idx);
    data_diameters_sf{idx, 6} = vmd90_sf(idx);
    data_diameters_sf{idx, 7} = ch_sf(idx);
    data_diameters_sf{idx, 8} = mean_diameter_sf(idx);
    data_diameters_sf{idx, 9} = std_diameter_sf(idx);
end

close all
return
