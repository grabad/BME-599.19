%% amri_fmri_smooth 
% spatial smoothing using a 3-D guassian kernel
%
% Usage
%   output = amri_fmri_smooth(input,dxyz,fwhm)
% Inputs 
%    img: a volume or time series of volumes
%   dxyz: voxel size (in mm)
%   fwhm: full-width-half-maximum (in mm)
% Outputs 
%   oimg: smoothed images
%
% Version
%    0.04
%
% History
%    0.00 - ZMLIU - 06/06/2013 - Original version calling spm_smooth
%    0.01 - ZMLIU - 06/12/2014 - use a fft-based convolution
%    0.02 - ZMLIU - 09/22/2014 - do not round when defining gaussian kernel
%    0.03 - ZMLIU - 03/18/2015 - fix a bug related to mask
%    0.04 - ZMLIU - 10/01/2023 - convolution with a 3D gaussian kernel

%%
function oimg = amri_fmri_smooth(img,dxyz,fwhm,varargin)

if nargin<1
    help('amri_fmri_smooth');
    return
end

if isnan(fwhm) || fwhm==0
    oimg=img;
    return;
end

if isscalar(fwhm)
    fwhm=ones(3,1)*fwhm;
end

mask=nan; 

for i = 1:2:size(varargin,2) 
    Keyword = varargin{i};
    Value   = varargin{i+1};
    if ~ischar(Keyword)
        fprintf('amri_fmri_smooth(): keywords must be strings\n'); 
        continue;
    end
    if strcmpi(Keyword,'mask')
        mask=double(Value);
    end
end

[nx,ny,nz,nt]=size(img);
if ~isnan(mask) 
    if size(mask,1)~=nx || ...
       size(mask,2)~=ny || ...
       size(mask,3)~=nz
        fprintf('amri_fmri_smooth(): dimension mismatch\n');
        fprintf('amri_fmri_smooth(): ignore the input mask\n');
        return;
    else
        mask=std(img,0,4)>0;
    end
else
    mask=std(img,0,4)>0;
end

fwhm=double(fwhm(:)');
dxyz=double(dxyz(:)');

% dx=dxyz(1);
% dy=dxyz(2);
% dz=dxyz(3);

% Convert FWHM to sigma (standard deviation)
sigma = fwhm / (2 * sqrt(2 * log(2)));
sigma_voxels = sigma ./ dxyz; 

% Calculate kernel size in voxel units
k_size_voxels = ceil(3 * sigma_voxels);

[x, y, z] = ndgrid(-k_size_voxels(1):k_size_voxels(1), ...
                   -k_size_voxels(2):k_size_voxels(2), ...
                   -k_size_voxels(3):k_size_voxels(3));

% Define Gaussian kernel for x, y, and z in voxel units
gauss_kernel_x = exp(-x.^2 / (2*sigma_voxels(1)^2));
gauss_kernel_y = exp(-y.^2 / (2*sigma_voxels(2)^2));
gauss_kernel_z = exp(-z.^2 / (2*sigma_voxels(3)^2));

% 3D Gaussian kernel in voxel units
gaussian_kernel_voxels = gauss_kernel_x .* gauss_kernel_y .* gauss_kernel_z;

% Normalize the kernel
gaussian_kernel_voxels = gaussian_kernel_voxels / sum(gaussian_kernel_voxels(:));

% Apply the Gaussian kernel at each time point
oimg = zeros(size(img));
for t = 1:nt
    % Mask the current time frame
    img_frame = img(:,:,:,t);
    img_frame(mask == 0) = 0;
    % Convolve with the Gaussian kernel
    convolved_frame = convn(img_frame, gaussian_kernel_voxels, 'same');
    % Replace voxels outside the mask with original values
    convolved_frame(mask == 0) = img_frame(mask == 0);    
    % Assign to the output image
    oimg(:,:,:,t) = convolved_frame;    
end




