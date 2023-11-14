clear
close all

load("data.mat")
window = 0.0005;
channel = 5;

data_length = length(data);
time_length = data_length/fs;
time_span = 0:1/fs:time_length-(1/fs);

%%
means = mean(data, 2);
stdevs = std(data, 0, 2);

peaks = find_spike_peaks(data(channel,:), 50, 75, fs, window);

figure()
hold on
plot(time_span, data(channel,:))
scatter(time_span(peaks), data(channel, peaks), 10, 'filled')
xlim([0 time_span(end)])
xlabel('Time (s)')
ylabel('Voltage (uV)')
title(['Neural Spike Data w/ Peaks, Channel ', num2str(channel)])

figure()
hold on
plot(time_span, data(channel,:))
scatter(time_span(peaks), data(channel, peaks), 10, 'filled')
xlim([0 100])
xlabel('Time (s)')
ylabel('Voltage (uV)')
title(['Neural Spike Data w/ Peaks, Channel ', num2str(channel), ' (First 100s)'])