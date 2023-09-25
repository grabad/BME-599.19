clear
close all

load("proj01_data.mat")
timespan = 0 : 1/srate : size(meg_data,1)/srate - 1/srate;


filtered_meg_data = zeros(size(meg_data));

for i=1:14
    filtered_data = amri_sig_filtfft(meg_data(:, i), srate, 0.5, 59);
    filtered_meg_data(:,i) = filtered_data;
end

writematrix(filtered_meg_data, 'filtered_meg_data.csv')