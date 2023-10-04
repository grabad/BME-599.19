%% load data and retrive information

% load the data file func.nii (Nifti format)
nii = amri_file_loadnii('func.nii');

% visualize the data using amri_fmri_sliceview.m 
% try help amri_fmri_sliceview for instructions about this tool
% explore this tool, e.g. clicking on a location to see the time series
amri_fmri_sliceview(nii.img);

% nii includes hdr (header) and img (4-D array)

% nx: the x dimension
% ny: the y dimension
% nz: the z dimension
% nt: the t dimension 
[nx,ny,nz,nt] = size(nii.img);

% the header includes the voxel size and TR (sampling interval)
dx = nii.hdr.pixdim(2); % the unit is 1/10 mm
dy = nii.hdr.pixdim(3); % the unit is 1/10 mm
dz = nii.hdr.pixdim(4); % the unit is 1/10 mm
dt = nii.hdr.pixdim(5); % the unit is seconds

%% spatial smoothing by convolving a 3-D gaussian kernel with a fwhm=10 
% (the unit is 1/10mm)
fwhm=10; 

% nii.img = amri_fmri_smooth(nii.img, [dx dy dz], fwhm); 

%% reshape nii.img as img, which is a 2-D array 
% the number of rows is the number of spatial locations 
% the number of columns is the number of time samples

img = reshape(nii.img,nx*ny*nz,nt);

%% run svd with img

% [U,S,V] = svd(img,'econ');


%% for each location, detrend the signal by removing the 3rd polynomial fit 

for i = 1 : size(img,1)
%    img(i,:) = amri_sig_detrend(img(i,:), 3, dt);
end

%% run svd with img again 


%% 
load stim.mat;
load hrf.mat
