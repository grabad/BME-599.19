function spike = spike_shape(data, spike_index, fs, window)
    start = round(spike_index - fs*window);
    if start < 1
        start = 1;
    end

    stop = round(spike_index + fs*window);
    if stop > length(data)
        stop = length(data);
    end

    spike = data(start:stop);
end