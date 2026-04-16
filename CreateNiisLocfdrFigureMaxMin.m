% =============================================================================
% CreateNiisLocfdrFigureMaxMin.m
% =============================================================================
% Part 2 of 3: locfdr + CCA weights → NIfTI brain maps → BrainNet renders
%
% Purpose:
%   Read the MATLAB export from ExportLocfdrCcaWeights.R containing locfdr-
%   thresholded CCA weights (3244 voxels × 49 frequency bins, 0.78–19.1 Hz).
%   For each of four frequency bands (Delta, Theta, Alpha, Beta):
%     1. Create signed map: strongest significant weight per voxel
%     2. Create max map: most positive significant weight per voxel
%     3. Create min map: most negative significant weight per voxel
%   Export all three as NIfTI files for BrainNet visualization.
%
% Input:
%   cca_wx_locfdr_component0_1.mat -- output from ExportLocfdrCcaWeights.R
%
% Output:
%   CCA_locfdr_figure_maxmin_maps0_1/
%     CCA_locfdr_fdrN_Delta.nii, .._max.nii, .._min.nii
%     CCA_locfdr_fdrN_Theta.nii, .._max.nii, .._min.nii
%     CCA_locfdr_fdrN_Alpha.nii, .._max.nii, .._min.nii
%     CCA_locfdr_fdrN_Beta.nii, .._max.nii, .._min.nii
%
% Workflow:
%   1. ExportLocfdrCcaWeights.R → cca_wx_locfdr_component0_1.mat
%   2. CreateNiisLocfdrFigureMaxMin.m (this file) → NIfTI files
%   3. PlotBrainNetLocfdrDisplayBatch4ViewsBipolar.m → PNG renders
%
% Dependencies:
%   - J2img3d_nors_sg3244.m (NIfTI writer for 3244-voxel source space)
%
% =============================================================================

clear; clc;

bands = [1.5 3.5; 3.8 7.5; 7.9 12.5; 12.9 19.1];
bnames = {'Delta','Theta','Alpha','Beta'};
f0 = 0.78;
fres = 0.3906;
scale_factor = 1e4;

pathd = 'D:\OneDrive - CCLAB\Maria\Mediatio NeuroEpo\NeuroEPOMDS-main\To upload on CC Lab Github\';
mat_fname = fullfile(pathd, 'cca_wx_locfdr_component0_1.mat');
out_dir = fullfile(pathd, 'CCA_locfdr_figure_maxmin_maps0_1');

if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

S = load(mat_fname);

required_fields = {'wx_sig_weights', 'fdr_threshold', 'n_voxels', 'n_freqs'};
for k = 1:numel(required_fields)
    if ~isfield(S, required_fields{k})
        error('Missing field "%s" in %s.', required_fields{k}, mat_fname);
    end
end

n_voxels = double(S.n_voxels);
n_freqs = double(S.n_freqs);

if numel(S.wx_sig_weights) ~= n_voxels * n_freqs
    error('wx_sig_weights does not match n_voxels x n_freqs.');
end

Wsig = reshape(double(S.wx_sig_weights), n_voxels, n_freqs);
freq = f0 + (0:n_freqs-1) * fres;

bind = round((bands - f0) ./ fres) + 1;
bind(:,1) = max(bind(:,1), 1);
bind(:,2) = min(bind(:,2), n_freqs);

fdr_pct = round(double(S.fdr_threshold) * 100);

for k = 1:size(bind,1)
    idx = bind(k,1):bind(k,2);
    Wband = Wsig(:, idx);

    % Preserve the existing signed map: one value per voxel from the
    % strongest significant bin in absolute value within the band.
    [~, peak_idx] = max(abs(Wband), [], 2);
    linear_idx = sub2ind(size(Wband), (1:size(Wband,1))', peak_idx);
    Jsigned = Wband(linear_idx);
    Jsigned(all(Wband == 0, 2)) = 0;

    % Match the max/min spectrum figure at the band level by selecting the
    % single significant frequency with the largest positive max and the
    % single significant frequency with the most negative min.
    Wmax_sig = max(Wband, [], 1);
    Wmin_sig = min(Wband, [], 1);

    no_pos_col = all(Wband <= 0, 1);
    no_neg_col = all(Wband >= 0, 1);
    Wmax_sig(no_pos_col | Wmax_sig == 0) = NaN;
    Wmin_sig(no_neg_col | Wmin_sig == 0) = NaN;

    Jmax = zeros(n_voxels, 1);
    Jmin = zeros(n_voxels, 1);

    if any(~isnan(Wmax_sig))
        [~, rel_max_idx] = max(Wmax_sig);
        max_col = Wband(:, rel_max_idx);
        Jmax(max_col > 0) = max_col(max_col > 0);
        max_freq_hz = freq(idx(rel_max_idx));
    else
        max_freq_hz = NaN;
    end

    if any(~isnan(Wmin_sig))
        [~, rel_min_idx] = min(Wmin_sig);
        min_col = Wband(:, rel_min_idx);
        Jmin(min_col < 0) = min_col(min_col < 0);
        min_freq_hz = freq(idx(rel_min_idx));
    else
        min_freq_hz = NaN;
    end

    signed_nii = fullfile(out_dir, sprintf('CCA_locfdr_fdr%d_%s.nii', fdr_pct, bnames{k}));
    max_nii = fullfile(out_dir, sprintf('CCA_locfdr_fdr%d_%s_max.nii', fdr_pct, bnames{k}));
    min_nii = fullfile(out_dir, sprintf('CCA_locfdr_fdr%d_%s_min.nii', fdr_pct, bnames{k}));

    J2img3d_nors_sg3244(Jsigned * scale_factor, signed_nii, 0);
    J2img3d_nors_sg3244(Jmax * scale_factor, max_nii, 0);
    J2img3d_nors_sg3244(Jmin * scale_factor, min_nii, 0);

    fprintf('%s band\n', bnames{k});
    fprintf('  Signed map voxels: %d\n', nnz(Jsigned ~= 0));
    fprintf('  Max map voxels: %d', nnz(Jmax > 0));
    if isnan(max_freq_hz)
        fprintf(' (no positive significant frequency in band)\n');
    else
        fprintf(' at %.4f Hz\n', max_freq_hz);
    end
    fprintf('  Min map voxels: %d', nnz(Jmin < 0));
    if isnan(min_freq_hz)
        fprintf(' (no negative significant frequency in band)\n');
    else
        fprintf(' at %.4f Hz\n', min_freq_hz);
    end
end

fprintf('Saved figure-style locfdr NIfTI files in %s\n', out_dir);