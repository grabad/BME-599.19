clear
close all

tic
%% Q1
% Load the image file. Specify the image size in terms of the number of spatial locations and time points. 
% What is the size of a voxel? What is the sampling interval? 
% How long is the fMRI recording in terms of seconds and minutes? 
% (Hints: the size of the voxel is specified in nii.hdr.pixdim(2:4) with a unit of 1/10 mm; the time interval is 
% specified in nii.hdr.pixdim(5) with a unit of seconds. To display the data, using amri_fmri_sliceview.m). 

img = amri_file_loadnii('func.nii');
hdr = img.hdr;
data = img.img;

%amri_fmri_sliceview(data)

%% Q2
% Spatially filter or smooth the images by applying a 3-D Gaussian filter with a full-width-at-half-maximum (FWHM=1mm or 10 * 1/10mm). 
% This can be done either by using 3-D convolution or using pairs of 3-D fourier transform and inverse fourier transform. 
% Use the smoothed data for subsequent analyses. 

filt_data = amri_fmri_smooth(data, hdr.pixdim(2)/10, 10/10);
%amri_fmri_sliceview(filt_data)

%% Q3
% Run SVD with the dataset. Plot the singular values. Report and discuss your observations.
data_size = size(data);
filt_data_flat = reshape(filt_data, [prod(data_size(1:3)), 590]);
[U, S, V] = svd(filt_data_flat, 'econ');

figure()
svs = diag(S);
plot(svs)

title('Singular Values: Filtered Data')
xlabel('Singular Value Number')
ylabel('Singular Value Magnitude')
xlim([0 100])

%% Q4
% For each voxel, remove a slow drift or detrend by regressing out a 3rd-order polynomial function that fits the slow detrend. 
% Use the detrended data for the subsequent analyses. 

detrend_data = detrend(filt_data_flat', 3)';

%% Q5
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
subplot(2,1,1)
hold on
plot(Ud(:, 1:3))

title('Left Singular Vectors: Detrended Data')
legend('LSV 1', 'LSV 2', 'LSV 3')
xticklabels([])

subplot(2,1,2)
hold on
plot(Vd(:, 1:3))

title('Right Singular Vectors: Detrended Data')
legend('RSV 1', 'RSV 2', 'RSV 3', Location='southeast')
xticklabels([])

%% Q6
% Keep the largest K principal components such that these components explain the 90% of the total variance and discard other components 
% from the data. What is the value of K? 
warning('off')
[coeff, score, latent, tsquared, explained] = pca(detrend_data);

sum_var = 0;
for i=1:590
    if sum_var < 90
        sum_var = sum_var + explained(i);
        K = i;
    else
        break
    end
end

%% Q7
% For each voxel, standardize the signal, such that the mean is zero and standard deviation is 1. 
standard_data = zeros(size(detrend_data));

for i=1:size(detrend_data, 1)
    row_data = detrend_data(i,:);
    standard_data(i,:) = (row_data -  mean(row_data))./std(row_data);
end

standard_data(isnan(standard_data)) = 0;

%% Q8
% Load stim.mat. Two types of stimuli are stored as the two rows in stim_block. Plot each stimulus as a time series. 
% The stimulus time series is sampled in synchronization with the fMRI. 
load('stim.mat');
time = linspace(0,708,590);

figure()
subplot(2,1,1)
plot(time, stim_block(1,:))

title('Stimulus #1')
xlabel('Time [sec]')
ylabel('Stimulus Magnitude [au]')
xlim([0 708])

subplot(2,1,2)
plot(time, stim_block(2,:))

title('Stimulus #2')
xlabel('Time [sec]')
ylabel('Stimulus Magnitude [au]')
xlim([0 708])

%% Q9
% Convolve each of the two stimuli with a hemodynamic response function (HRF), which can be loaded from hrf.mat in the folder. 
% Demean each Consider and plot the results as two regressors for subsequent analyses.
load('hrf.mat')
stim1 = conv(stim_block(1,:), hrf, 'same');
stim2 = conv(stim_block(2,:), hrf, 'same');

reg1 = detrend(stim1,'constant');
reg2 = detrend(stim2,'constant');

figure()
subplot(2,1,1)
plot(reg1)

title('Regressor #1')
ylabel('Regressor Magnitude [au]')
xlim([0 590])

subplot(2,1,2)
plot(reg2)

title('Regressor #2')
ylabel('Regressor Magnitude [au]')
xlim([0 590])

%% Q10
% Describe a linear regression model that predicts each voxel’s signal as a linear combination of the two regressors obtained in 9). 
X = [reg1', reg2'];
betas = zeros([2, 590]);

for i=1:size(standard_data, 1)
    Y = standard_data(i,:)';
    betas(:, i) = (X'*X)\X'*Y;
end

%% Q11/12
% For each voxel, evaluate the model performance using F statistic. 
% For each voxel, evaluate the significance that each regressor explains the signal.

fStats = zeros([124800,1]);
pVals = zeros([124800,1]);

for i = 1:size(standard_data, 1)
    voxel = standard_data(i, :)';
    residuals = voxel - X*betas(:, i);

    SSR = sum((X*betas(:, i) - mean(voxel)).^2);
    SSE = sum(residuals.^2);
    
    % Calculate degrees of freedom
    df_SSR = 2 - 1;
    df_SSE = 590 - 2;
    
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

%% Q13
% Show a map of the significant voxels above a significance level of 0.05. 

pVal_map = reshape(pVals,[40,78,40]);
binary_map = pVal_map < 0.05;

[x, y, z] = meshgrid(1:78, 1:40, 1:78);
significant_voxels = find(binary_map);

figure()
scatter3(x(significant_voxels), y(significant_voxels), z(significant_voxels), 'filled');
title('Voxels Significant to 0.05')

%% Q14
% Show the map again after correction for multiple comparisons using false discovery rate < 0.05. 
corrected_pVals = fdr_bh(pVals);

binary_corrected_pVals_map = reshape(corrected_pVals, [40, 78, 40]);
significant_voxels_fdr = find(binary_corrected_pVals_map);

figure()
scatter3(x(significant_voxels_fdr), y(significant_voxels_fdr), z(significant_voxels_fdr), 'filled');
title({'Voxels Significant to 0.05', 'Corrected for Multiple Comparisons'})
toc