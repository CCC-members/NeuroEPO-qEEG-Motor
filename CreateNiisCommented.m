% CreateNiisCommented.m
% This script processes single-cca scores, applies frequency band selection, and converts the results to 3D NIfTI images.

% Define frequency bands in Hz
bands = [1.5 3.5; 3.8 7.5; 7.9 12.5; 12.9 19.1];
bnames = {'Delta', 'Theta', 'Alpha', 'Beta'}; % Corresponding band names
fres = 0.3906; % Frequency resolution
bind = ceil(bands ./ fres); % Convert bands to index ranges

% Define the base directory for loading and saving data
pathd = '';

% Load the single CCA scores dataset
load([pathd filesep 'single_cca_scores.mat']);


% Process the single CCA score dataset
% nii_fname = [pathd 'single_cca_scores2.nii'];
nii_fnameDelta = [pathd 'single_cca_scores2_Delta.nii'];
nii_fnameTheta = [pathd 'single_cca_scores2_Theta.nii'];
nii_fnameAlpha = [pathd 'single_cca_scores2_Alpha.nii'];
nii_fnameBeta = [pathd 'single_cca_scores2_Beta.nii'];
single_cca_score2 = reshape(single_cca_score2, 3244, 49) * 10000; % Reshape and scale

Jb = []; % Reset the storage matrix

for k=1:size(bind,1)
    % Compute log power for the maximum signal within each frequency band
    Jb(:,k) = log(max(single_cca_score2(:, bind(k,1):bind(k,2)), [], 2).^2);
    
    % Thresholding: set values below the 95th percentile to -120
    th1 = prctile(Jb(:,k),95);
    Jb(Jb(:,k) < th1,k) = -120;
end

% Convert the processed data into a 3D image
% J2img3d_nors_sg3244(Jb, nii_fname, 0);
J2img3d_nors_sg3244(Jb(:,1), nii_fnameDelta, 0);
J2img3d_nors_sg3244(Jb(:,2), nii_fnameTheta, 0);
J2img3d_nors_sg3244(Jb(:,3), nii_fnameAlpha, 0);
J2img3d_nors_sg3244(Jb(:,4), nii_fnameBeta, 0);

