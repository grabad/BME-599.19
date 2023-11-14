function peaks = find_spike_peaks(data, top_thresh, bottom_thresh, fs, window)
    [~, peaks_p] = findpeaks(data, 'MinPeakHeight', top_thresh, 'MinPeakDistance', fs*window, 'MinPeakProminence', 30);
    [~, peaks_n] = findpeaks(-data, 'MinPeakHeight', bottom_thresh, 'MinPeakDistance', fs*window, 'MinPeakProminence', 30);
    peaks = sort([peaks_p, peaks_n]);

    i = 1;
    while i < length(peaks)
        if peaks(i+1) - peaks(i) < fs*window
            if abs(data(peaks(i))) > abs(data(peaks(i+1)))
                peaks(i+1) = [];
            else
                peaks(i) = [];
            end
        else
        i = i + 1;
        end
    end
end