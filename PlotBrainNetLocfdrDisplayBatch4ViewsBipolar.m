% =============================================================================
% PlotBrainNetLocfdrDisplayBatch4ViewsBipolar.m
% =============================================================================
% Part 3 of 3: locfdr + CCA weights → NIfTI brain maps → BrainNet renders
%
% Purpose:
%   Read the NIfTI files output by CreateNiisLocfdrFigureMaxMin.m and render
%   each as a 4-view BrainNet figure (left sagittal, right sagittal,
%   axial top, axial bottom) using a symmetric blue-white-red-yellow
%   diverging colormap anchored at zero.
%   Filters to only the signed (bipolar) maps, skipping per-band max/min
%   specializations that were used for spectrum figures.
%
% Input:
%   CCA_locfdr_figure_maxmin_maps0_1/
%     CCA_locfdr_fdrN_Delta.nii
%     CCA_locfdr_fdrN_Theta.nii
%     CCA_locfdr_fdrN_Alpha.nii
%     CCA_locfdr_fdrN_Beta.nii
%
% Output:
%   BrainNet_Locfdr_Display_4Views_Output0_1_BipolarOnly/
%     [for each band]/
%       *_left_sagittal.png
%       *_right_sagittal.png
%       *_axial_top.png
%       *_axial_bottom.png
%
% Workflow:
%   1. ExportLocfdrCcaWeights.R → cca_wx_locfdr_component0_1.mat
%   2. CreateNiisLocfdrFigureMaxMin.m → NIfTI files
%   3. PlotBrainNetLocfdrDisplayBatch4ViewsBipolar.m (this file) → PNG renders
%
% Dependencies:
%   - BrainNet Viewer (https://www.nitrc.org/projects/bnv/)
%   - Custom helpers: reference_diverging_colormap, compute_symmetric_range,
%     adjust_brainnet_colormap, set_brainnet_rest_flag,
%     adjust_sagittal_colorbar_layout, export_brainnet_png
%
% =============================================================================

clear; clc; close all;

%% Paths
brainnet_dir = 'D:\OneDrive\Downloads\BrainNetViewer_20191031';
surf_file = 'D:\OneDrive\Downloads\BrainNetViewer_20191031\Data\SurfTemplate\BrainMesh_ICBM152.nv';

pathd = 'D:\OneDrive - CCLAB\Maria\Mediatio NeuroEpo\NeuroEPOMDS-main\To upload on CC Lab Github\';
input_dir = fullfile(pathd, 'CCA_locfdr_figure_maxmin_maps0_1');
cfg_template = fullfile(pathd, 'cfg_left_lateral.mat');
out_root = fullfile(pathd, 'BrainNet_Locfdr_Display_4Views_Output0_1_BipolarOnly');

if ~exist(brainnet_dir, 'dir')
    error('BrainNet directory not found: %s', brainnet_dir);
end

if ~exist(surf_file, 'file')
    error('Surface file not found: %s', surf_file);
end

if ~exist(input_dir, 'dir')
    error('Input directory not found: %s', input_dir);
end

if ~exist(cfg_template, 'file')
    error('Config template not found: %s', cfg_template);
end

if ~exist(out_root, 'dir')
    mkdir(out_root);
end

addpath(genpath(brainnet_dir));
addpath(pathd);

nii_files = dir(fullfile(input_dir, '*.nii'));
if isempty(nii_files)
    error('No NIfTI files found in %s', input_dir);
end

keep_idx = ~contains({nii_files.name}, '_max.nii') & ~contains({nii_files.name}, '_min.nii');
nii_files = nii_files(keep_idx);
if isempty(nii_files)
    error('No combined signed NIfTI files found in %s', input_dir);
end

view_specs = {
    'left_sagittal',  -90,    0;
    'right_sagittal',  90,    0;
    'axial_top',        0,   90;
    'axial_bottom',  -180,  -90;
};

base_cfg = load(cfg_template, 'EC');
reference_map = reference_diverging_colormap(1001);

for file_idx = 1:numel(nii_files)
    nii_name = nii_files(file_idx).name;
    nii_file = fullfile(input_dir, nii_name);
    [~, nii_stem] = fileparts(nii_name);
    map_outdir = fullfile(out_root, nii_stem);

    if ~exist(map_outdir, 'dir')
        mkdir(map_outdir);
    end

    fprintf('Processing %s\n', nii_name);

    for view_idx = 1:size(view_specs, 1)
        view_name = view_specs{view_idx, 1};
        view_az = view_specs{view_idx, 2};
        view_el = view_specs{view_idx, 3};

        EC = base_cfg.EC;
        EC.lot.view = 1;
        EC.lot.view_direction = 4;
        EC.lot.view_az = view_az;
        EC.lot.view_el = view_el;
        EC.img.width = 2200;
        EC.img.height = 1800;
        EC.img.dpi = 300;

        [negative_max, negative_min, positive_min, positive_max] = compute_symmetric_range(nii_file);

        EC.vol.display = 1;
        EC.vol.null = [0.98 0.98 0.98];
        EC.vol.nx = negative_max;
        EC.vol.nn = negative_min;
        EC.vol.pn = positive_min;
        EC.vol.px = positive_max;
        EC.vol.color_map = 24;
        EC.vol.cmstring = 'custom_symmetric_reference_scale';
        EC.vol.CMt = reference_map;
        EC.vol.CM = adjust_brainnet_colormap(reference_map, EC.vol.null, EC.vol.nx, EC.vol.nn, EC.vol.pn, EC.vol.px);

        temp_cfg_file = fullfile(map_outdir, sprintf('%s_%s_temp_cfg.mat', nii_stem, view_name));
        out_png = fullfile(map_outdir, sprintf('%s_%s.png', nii_stem, view_name));

        save(temp_cfg_file, 'EC');

        fprintf('  Rendering %s...\n', view_name);
        set_brainnet_rest_flag();

        h = BrainNet_MapCfg(surf_file, nii_file, temp_cfg_file);
        drawnow;

        if contains(view_name, 'sagittal')
            adjust_sagittal_colorbar_layout(h);
            drawnow;
        end

        export_brainnet_png(h, out_png, EC.img.dpi);

        if ishghandle(h)
            close(h);
        end

        if ~exist(out_png, 'file')
            error('BrainNet did not create the expected output file: %s', out_png);
        end

        delete(temp_cfg_file);
        fprintf('  Saved %s\n', out_png);
    end
end

fprintf('Finished rendering all bipolar maps to %s\n', out_root);

function new_colormap = adjust_brainnet_colormap(original_colormap, null_color, negative_max, negative_min, positive_min, positive_max)
temp_colormap = original_colormap;
sample_idx = floor(linspace(1, size(temp_colormap, 1) + 0.9999, 1000));
original_colormap = temp_colormap(sample_idx, :);
new_colormap = repmat(null_color, [1000 1]);

negative_color_segment = fix(1000 * (negative_min - negative_max) / (positive_max - negative_max));
if negative_color_segment == 0
    negative_color_segment = 1;
end

negative_index = round(linspace(1, 500, negative_color_segment));
new_colormap(1:negative_color_segment, :) = original_colormap(negative_index, :);

positive_color_segment = fix(1000 * (positive_max - positive_min) / (positive_max - negative_max));
if positive_color_segment == 0
    positive_color_segment = 1;
end

positive_index = round(linspace(501, 1000, positive_color_segment));
new_colormap(end - positive_color_segment + 1:end, :) = original_colormap(positive_index, :);
end

function set_brainnet_rest_flag()
evalin('base', ['global FLAG; ' ...
    'if ~exist(''FLAG'',''var'') || ~isstruct(FLAG), FLAG = struct(); end; ' ...
    'FLAG.IsCalledByREST = 1;']);
end

function adjust_sagittal_colorbar_layout(fig_handle)
colorbars = findall(fig_handle, 'Type', 'ColorBar');
if isempty(colorbars)
    colorbars = findall(fig_handle, 'Tag', 'Colorbar');
end

if isempty(colorbars)
    return;
end

main_axes = findall(fig_handle, 'Type', 'Axes');

axes_tags = get(main_axes, 'Tag');
if ischar(axes_tags)
    axes_tags = {axes_tags};
end

keep_axes = true(size(main_axes));
for k = 1:numel(main_axes)
    pos = main_axes(k).Position;
    is_colorbar_axis = strcmpi(axes_tags{k}, 'Colorbar');
    is_full_figure_axis = pos(1) <= 0.01 && pos(2) <= 0.01 && pos(3) >= 0.99 && pos(4) >= 0.99;
    keep_axes(k) = ~is_colorbar_axis && ~is_full_figure_axis;
end

main_axes = main_axes(keep_axes);

reserved_right_edge = 0.83;
colorbar_position = [0.90 0.10 0.03 0.30];

for k = 1:numel(main_axes)
    pos = main_axes(k).Position;
    right_edge = pos(1) + pos(3);

    if right_edge > reserved_right_edge
        pos(3) = max(0.01, reserved_right_edge - pos(1));
        main_axes(k).Position = pos;
    end
end

for k = 1:numel(colorbars)
    colorbars(k).Units = 'normalized';
    colorbars(k).Position = colorbar_position;
end
end

function export_brainnet_png(fig_handle, out_png, dpi)
paper_width = 2200 / dpi;
paper_height = 1800 / dpi;

set(fig_handle, 'Units', 'pixels');
set(fig_handle, 'Position', [100 100 2200 1800]);
set(fig_handle, 'PaperPositionMode', 'manual');
set(fig_handle, 'PaperUnits', 'inch');
set(fig_handle, 'PaperSize', [paper_width paper_height]);
set(fig_handle, 'PaperPosition', [0 0 paper_width paper_height]);
print(fig_handle, out_png, '-dpng', ['-r', num2str(dpi)]);
end

function cm = reference_diverging_colormap(m)
anchor_values = [-4.0 -3.5 -3.0 -2.5 -2.0 -1.0 -0.5 0.0 0.5 1.0 2.0 2.5 3.0 3.5 4.0];
anchor_colors = [
    0.0000 0.0500 0.4500
    0.0500 0.1800 0.7000
    0.0500 0.3000 0.8500
    0.0500 0.4500 0.9500
    0.0500 0.7000 0.9500
    0.4500 0.9000 1.0000
    0.8200 0.9700 0.9900
    0.9800 0.9800 0.9800
    1.0000 0.9400 0.9000
    1.0000 0.7800 0.6000
    1.0000 0.4200 0.2200
    1.0000 0.2800 0.1000
    0.9300 0.0800 0.0800
    0.9800 0.3300 0.0500
    1.0000 0.9300 0.0000
];

target_values = linspace(-4, 4, m);
cm = interp1(anchor_values, anchor_colors, target_values, 'pchip');
end

function [negative_max, negative_min, positive_min, positive_max] = compute_symmetric_range(nii_file)
vol_hdr = BrainNet_spm_vol(nii_file);
vol_img = BrainNet_spm_read_vols(vol_hdr);
nonzero_vals = vol_img(~isnan(vol_img) & vol_img ~= 0);

if isempty(nonzero_vals)
    positive_max = 1;
    positive_min = 0.02;
    negative_min = -0.02;
    negative_max = -1;
    return;
end

abs_limit = max(abs(nonzero_vals));

inner_fraction = 0.02;
inner_limit = max(abs_limit * inner_fraction, eps);

positive_min = inner_limit;
negative_min = -inner_limit;
positive_max = abs_limit;
negative_max = -abs_limit;
end