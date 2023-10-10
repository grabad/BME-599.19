clear
close all

% Load the image file. Specify the image size in terms of the number of spatial locations and time points. 
% What is the size of a voxel? What is the sampling interval? 
% How long is the fMRI recording in terms of seconds and minutes? 
% (Hints: the size of the voxel is specified in nii.hdr.pixdim(2:4) with a unit of 1/10 mm; the time interval is 
% specified in nii.hdr.pixdim(5) with a unit of seconds. To display the data, using amri_fmri_sliceview.m). 
img = amri_file_loadnii('func.nii');
hdr = img.hdr;
data = img.img;

%amri_fmri_sliceview(data)

% Spatially filter or smooth the images by applying a 3-D Gaussian filter with a full-width-at-half-maximum (FWHM=1mm or 10 * 1/10mm). 
% This can be done either by using 3-D convolution or using pairs of 3-D fourier transform and inverse fourier transform. 
% Use the smoothed data for subsequent analyses. 
filt_data_test = zeros(size(data));
sigma = 10/(2*sqrt(2*log(2)));

for i=1:size(data, 4)
    filt_data_test(:, :, :, i) = imgaussfilt3(data(:, :, :, i), sigma);
end

filt_data = amri_fmri_smooth(data, 0.5, 1);

%amri_fmri_sliceview(filt_data)

% Run SVD with the dataset. Plot the singular values. Report and discuss your observations.

filt_data_flat = reshape(filt_data, [124800 590]);
[U, S, V] = svd(filt_data_flat, 'econ');

figure()
svs = diag(S);
plot(svs)

title('Singular Values: Filtered Data')
xlabel('Singular Value Number')
ylabel('Singular Value Magnitude')
xlim([0 100])

% For each voxel, remove a slow drift or detrend by regressing out a 3rd-order polynomial function that fits the slow detrend. 
% Use the detrended data for the subsequent analyses. 

detrend_data = detrend(filt_data_flat', 3)';

% Run SVD with the detrended data. Plot the singular values. Visualize the first three left and right singular vectors. 
% (Hint: if the left singular vectors are in U and the right singular vectors are in V, where SVD(A) = U * S * V’. 
% Then the left singular vectors are spatial patterns. The right singular vectors are temporal patterns.)

[Ud, Sd, Vd] = svd(detrend_data, 'econ');

figure()
svsd = diag(Sd);
plot(svsd)

title('Singular Values: Detrended Data')
xlabel('Singular Value Number')
ylabel('Singular Value Magnitude')
xlim([0 100])

figure()
hold on
plot(Ud())


% Keep the largest K principal components such that these components explain the 90% of the total variance and discard other components 
% from the data. What is the value of K? 
[coeff, score, latent] = pca(detrend_data);
K_initial = size(latent);
intitial_sum = sum(latent);
ninety_percent = .9*intitial_sum;

sum_var = 0;
K = 0;
for i=1:K_initial
    if sum_var < ninety_percent
        sum_var = sum_var + latent(i);
        K = K + 1;
    else
        break
    end
end


% For each voxel, standardize the signal, such that the mean is zero and standard deviation is 1. 
standard_data = (detrend(detrend_data,'constant'))./std(detrend_data);

%%
% Load stim.mat. Two types of stimuli are stored as the two rows in stim_block. Plot each stimulus as a time series. 
% The stimulus time series is sampled in synchronization with the fMRI. 
load('stim.mat');
time = linspace(0,708,590);

figure()
plot(time, stim_block(1,:))

title('Stimuli #1')
xlabel('Time [sec]')
ylabel('Stimuli [unspecified]')
xlim([0 708])

figure()
plot(time, stim_block(2,:))

title('Stimuli #2')
xlabel('Time [sec]')
ylabel('Stimuli [unspecified]')
xlim([0 708])


% Convolve each of the two stimuli with a hemodynamic response function (HRF), which can be loaded from hrf.mat in the folder. 
% Demean each Consider and plot the results as two regressors for subsequent analyses.
load('hrf.mat')
stim1 = conv(stim_block(1,:),hrf);
stim2 = conv(stim_block(2,:),hrf);

demean_stim1 = detrend(stim1,'constant');
demean_stim2 = detrend(stim2,'constant');

figure()
plot(demean_stim1)
hold on
plot(demean_stim2)


% Describe a linear regression model that predicts each voxel’s signal as a linear combination of the two regressors obtained in 9). 

X = [demean_stim1(1,1:590); demean_stim2(1,1:590); ones([1,590])]';
Y = standard_data';

B = inv(X'*X)*X'*Y;

% For each voxel, evaluate the model performance using F statistic. 
% For each voxel, evaluate the significance that each regressor explains the signal.

fStats = zeros([124800,1]);
pVals = zeros([124800,1]);

for i = 1:124800
    voxel = Y(:,i);
    residuals = voxel - X*B(:,i);

    SSR = sum((X*B(:,i) - mean(voxel)).^2);
    SSE = sum(residuals.^2);
    
    % Calculate degrees of freedom
    df_SSR = 2;
    df_SSE = 590 - 3;
    
    % Calculate mean squares
    MSR = SSR / df_SSR;
    MSE = SSE / df_SSE;
    
    % Calculate the F-statistic
    fStat = MSR / MSE;
    fStats(i,1) = fStat;

    % Calculate the p-value
    pVal = 1 - fcdf(fStat, df_SSR, df_SSE);
    pVals(i,1) = pVal;
end

figure()
plot(pVals)

% Show a map of the significant voxels above a significance level of 0.05. 

pVal_map = reshape(pVals,[40,78,40]);
binary_map = pVal_map > 0.05;           % note the only way you can see anything on the 3d plot is if you look at large p vals not small

[x, y, z] = meshgrid(1:78, 1:40, 1:78);
significant_voxels = find(binary_map);

figure()
scatter3(x(significant_voxels), y(significant_voxels), z(significant_voxels), 'filled', 'MarkerFaceColor', 'b');


% Show the map again after correction for multiple comparisons using false discovery rate < 0.05. 

[pFDR, q] = mafdr(pVals);
pFDR_map = reshape(pFDR, [40, 78, 40]);

binary_fdr_map = pFDR_map < 0.05;
significant_voxels_fdr = find(binary_fdr_map);

figure()
scatter3(x(significant_voxels_fdr), y(significant_voxels_fdr), z(significant_voxels_fdr), 'filled', 'MarkerFaceColor', 'r');

