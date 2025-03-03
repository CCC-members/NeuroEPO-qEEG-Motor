function imgfull = J2img3d_nors_sg3244Commented(J, nii_fname, change_to_aubert)
% J2img3d_nors_sg3244: Converts a data matrix (J) into a 3D image and saves it as a NIfTI file if required.
%
% INPUTS:
%   J - The input matrix containing the data.
%   nii_fname - The name of the output NIfTI file (optional).
%   change_to_aubert - A flag to apply a coordinate transformation (optional).
%
% OUTPUT:
%   imgfull - The generated 3D image.

% Load the NIfTI file containing the reference brain image
nii = load_nii('D:\User1 OneDrive\OneDrive - CCLAB\Maria\Mediatio NeuroEpo\NeuroEPOMDS-main\To upload on CC Lab Github\aal.nii');

% Load the coordinate mapping file
xyz = load('D:\User1 OneDrive\OneDrive - CCLAB\Maria\Mediatio NeuroEpo\NeuroEPOMDS-main\To upload on CC Lab Github\Nors_sg3244.xyz');
xyz = round(xyz); % Round coordinates to the nearest integer

% Adjust coordinate system
xyz = xyz(: , [1 3 2]); % Swap axes for compatibility
xyz(:, 2) = size(nii.img,2) - xyz(:, 2); % Invert y-coordinates
xyz(:, 1) = size(nii.img,1) - xyz(:, 1); % Invert x-coordinates

% Initialize an empty 4D image array
imgfull = zeros([size(nii.img) size(J,2)]);

% Replace NaN values in J with zeros
ii = find(isnan(J));
J(ii) = 0;

% Create a 3D grid for interpolation
[xq,yq,zq] = meshgrid(1:size(nii.img,2), 1:size(nii.img,1), 1:size(nii.img,3));

% Process each column in J separately
for k = 1:size(J,2)
    fprintf('.'); % Print progress indicator
    
    y = J(:,k);
    ii = find(y); % Find nonzero elements
    
    if length(ii) > 4
        img = zeros(size(nii.img)); % Initialize empty image
        
        % Assign values to the corresponding coordinates
        for h=1:length(ii)
            img(xyz(ii(h),1),xyz(ii(h),2),xyz(ii(h),3)) = y(ii(h));
        end
        
        % Perform interpolation on nonzero values
        ii = find(img);
        [xx,yy,zz] = ind2sub(size(img), ii);
        v = img(ii);
        img = griddata(yy,xx,zz,v,xq,yq,zq);
        imgfull(:,:,:,k) = img;
    end
end
fprintf('\n'); % Move to the next line

% Save the generated image as a NIfTI file if a filename is provided
if nargin > 1
    nii.img = imgfull;
    nii.hdr.dime.datatype = 16; % Set data type to float
    nii.hdr.dime.bitpix = 16; % Set bit depth
    save_nii(nii, nii_fname);
end

% Apply coordinate transformation if requested
if (nargin>2) && change_to_aubert
    xyz = xyz(: , [1 3 2]); % Swap axes again
    xyz(:, 2) = size(nii.img,2) - xyz(:, 2);
    xyz(:, 1) = size(nii.img,1) - xyz(:, 1);
end
