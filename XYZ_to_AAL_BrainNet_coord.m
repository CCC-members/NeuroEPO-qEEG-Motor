function XYZ_to_AAL_BrainNet_coord(x, save_file)
% x: matrix of nsources x nregs with the activations
% regs: index of the regs (colums) to be generated in the plot
% save_file: file where the images will be saved
%Author: Jorge Bosch Bayard
%% loading relevant data
addpath('')% for create_plot_4_BrainNet

fname_voxelscod = 'D:\FULEAH\MDS\AAL_NORSGC_index_3244.ind'; 
atlas_name = 'D:\FULEAH\MDS\aal.cod';
min_percent = 0; %porciento
ng = size(x,1);

addpath('D:\FULEAH\MDS\NIfTI_20140122')% add path for .nii reader 
%(Download at https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image)
nii = load_nii('D:\FULEAH\MDS\NITRC-multi-file-downloads\Data\ExampleFiles\AAL90\aal.nii');

xyz = load('D:\FULEAH\MDS\Nors_sg3244.xyz');
xyz = round(xyz);
%% Change XYZ coordinates reference system to AAL BrainNet
xyz = xyz(: , [1 3 2]);
xyz(:, 2) = size(nii.img,2) - xyz(:, 2);
xyz(:, 1) = size(nii.img,1) - xyz(:, 1);
full_img = zeros(size(nii.img));
for k=1:size(x,2) %The regions are interpolated one by one to avoid the interpolation joining them
    y = x(:,k);
    ii = find(y);
    img = zeros(size(nii.img));
    for h=1:length(ii)
        img(xyz(ii(h),1),xyz(ii(h),2),xyz(ii(h),3)) = y(ii(h));
    end
    ii = find(img);
    [xx,yy,zz] = ind2sub(size(img), ii);
    v = img(ii);
    [xq,yq,zq] = meshgrid(1:size(img,2), 1:size(img,1), 1:size(img,3));
    img = griddata(yy,xx,zz,v,xq,yq,zq);
    ii = find(~isnan(img));
    full_img(ii) = img(ii);
end
%% define new .nii file
nii.img = full_img;
nii.hdr.dime.datatype = 4;
nii.hdr.dime.bitpix = 4;
%% save new .nii file
save_nii(nii, save_file)
