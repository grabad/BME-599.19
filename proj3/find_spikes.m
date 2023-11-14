function spikes = find_spikes(data, peaks, fs, window)
    spikes = zeros([length(peaks), round(2*fs*window) + 1]);

    for i=1:length(peaks)
        spikes(i,:) = spike_shape(data, peaks(i), fs, window);
    end
end