import numpy as np
from scipy.io import loadmat 
from scipy.signal import detrend
from scipy.fft import fft, ifft
import seaborn as sns
from itertools import cycle 
import matplotlib.pyplot as plt


# Define PlotAll function
def plotAll(timespan, data, dataLabels, figTitle, figXLabel, figYLabel):
    fig = plt.figure(figsize=(12, 14))
    t = fig.add_gridspec(7, 2, wspace=0.05, hspace=0.05)

    colors = plt.cm.plasma(np.linspace(0, 1, 15))
    color_cycle = cycle(colors)

    timePad = 0.0042 * (timespan[-1] - timespan[0])

    for i in range(14):
        ax = fig.add_subplot(t[i])
        name = dataLabels[i]

        ax.plot(timespan, data[:, i], label=name, color=next(color_cycle))

        ax.set_xlim(timespan[0] - timePad, timespan[-1] + timePad)
        ax.set_ylim(np.min(data), np.max(data))

        if i != 12 and i != 13:
            ax.set_xticklabels([])

        if i % 2 == 1:
            ax.set_yticklabels([])

        ax.legend(loc='upper right')

    fig.text(0.5, 0.04, figXLabel, ha='center')
    fig.text(0.04, 0.5, figYLabel, va='center', rotation='vertical')
    fig.suptitle(figTitle, fontsize=16)

    plt.show()

# Define lowpass, bandpass, highpass filtering function
def fftFilt(data, lowCO=None, highCO=None):
    freqData = np.fft.fft(data)
    freqRange = np.fft.fftfreq(np.size(data), 1/300)

    if lowCO is not None:
        freqData[freqRange < lowCO] = 0
        
    if highCO is not None:
        freqData[freqRange > highCO] = 0
    
    filteredData = np.fft.ifft(freqData)
    
    return filteredData

def amri_sig_filtfft(ts, fs, lowcut, highcut, revfilt=0, trans=0.15):
    if not isinstance(ts, np.ndarray):
        raise ValueError('amri_sig_filtfft(): input data has to be a NumPy array')

    npts = len(ts)  # number of time points
    nfft = 2 ** int(np.ceil(np.log2(npts)))  # number of frequency points

    fv = np.linspace(0, fs / 2, nfft // 2 + 1)  # odd-sized frequency vector from 0 to Nyquist frequency
    fres = (fv[-1] - fv[0]) / (nfft // 2)  # frequency domain resolution

    # Remove the linear trend
    ts_old = ts.copy()
    ts = detrend(ts, type='linear')
    trend = ts_old - ts

    # Design frequency domain filter
    filter = np.ones(nfft)

    if (not np.isnan(lowcut) and lowcut > 0) and (np.isnan(highcut) or highcut <= 0):
        # Highpass filter
        idxl = int(round(lowcut / fres) + 1)
        idxlmt = int(round(lowcut * (1 - trans) / fres) + 1)
        idxlmt = max(idxlmt, 1)
        filter[0:idxlmt] = 0
        filter[idxlmt:idxl] = 0.5 * (1 + np.sin(-np.pi / 2 + np.linspace(0, np.pi, idxl - idxlmt + 1)))
        filter[nfft - idxl + 1:nfft] = filter[idxl::-1]

    elif (np.isnan(lowcut) or lowcut <= 0) and (not np.isnan(highcut) and highcut > 0):
        # Lowpass filter
        idxh = int(round(highcut / fres) + 1)
        idxhpt = int(round(highcut * (1 + trans) / fres) + 1)
        filter[idxh:idxhpt] = 0.5 * (1 + np.sin(np.pi / 2 + np.linspace(0, np.pi, idxhpt - idxh + 1)))
        filter[idxhpt:nfft // 2] = 0
        filter[nfft // 2 + 1:nfft - idxh + 1] = filter[nfft // 2:idxh - 1:-1]

    elif lowcut > 0 and highcut > 0 and highcut > lowcut:
        if revfilt == 0:  # Bandpass (revfilt == 0)
            transition = (highcut - lowcut) / 2 * trans
            idxl = int(round(lowcut / fres) + 1)
            idxlmt = int(round((lowcut - transition) / fres) + 1)
            idxh = int(round(highcut / fres) + 1)
            idxhpt = int(round((highcut + transition) / fres) + 1)
            idxl = max(idxlmt, 1)
            idxlmt = max(idxlmt, 1)
            idxh = min(nfft // 2, idxh)
            idxhpt = min(nfft // 2, idxhpt)
            filter[0:idxlmt] = 0
            filter[idxlmt:idxl] = 0.5 * (1 + np.sin(-np.pi / 2 + np.linspace(0, np.pi, idxl - idxlmt + 1)))
            filter[idxh:idxhpt] = 0.5 * (1 + np.sin(np.pi / 2 + np.linspace(0, np.pi, idxhpt - idxh + 1)))
            filter[idxhpt:nfft // 2] = 0
            filter[nfft - idxl + 1:nfft] = filter[idxl::-1]
            filter[nfft // 2 + 1:nfft - idxh + 1] = filter[nfft // 2:idxh - 1:-1]

        else:  # Bandstop (revfilt == 1)
            transition = (highcut - lowcut) / 2 * trans
            idxl = int(round(lowcut / fres) + 1)
            idxlmt = int(round((lowcut - transition) / fres) + 1)
            idxh = int(round(highcut / fres) + 1)
            idxhpt = int(round((highcut + transition) / fres) + 1)
            idxlmt = max(idxlmt, 1)
            idxlmt = max(idxlmt, 1)
            idxh = min(nfft // 2, idxh)
            idxhpt = min(nfft // 2, idxhpt)
            filter[idxlmt:idxl] = 0.5 * (1 + np.sin(np.pi / 2 + np.linspace(0, np.pi, idxl - idxlmt + 1)))
            filter[idxl:idxh] = 0
            filter[idxh:idxhpt] = 0.5 * (1 + np.sin(-np.pi / 2 + np.linspace(0, np.pi, idxhpt - idxh + 1)))
            filter[nfft - idxhpt + 1:nfft - idxlmt + 1] = filter[idxhpt::-1]

    else:
        raise ValueError('amri_sig_filtfft(): error in lowcut and highcut setting')

    X = fft(ts, nfft)  # FFT
    ts_new = np.real(ifft(X * filter, nfft))  # IFFT
    ts_new = ts_new[:npts]  # Truncate

    # Add back the linear trend
    ts_new = ts_new + trend

    return ts_new

# Import data and find time length in seconds and minutes
output = loadmat('proj1\proj01_data.mat')
meg_data = output['meg_data']
# meg_label from MATLAB array was coming up weird because it was a cell array, so we just rebuilt it
meg_label = ["MLF11", "MZF01", "MRF11", "MLC11", "MZC01", "MRC11", "MLP11", "MZP01", "MRP11", "MLO11", "MZO01", "MRO11", "MLT11", "MRT11"]
meg_locs = output['meg_locs']
srate = output['srate']

num_samples = np.shape(meg_data)[0]
srate_sec = 1/srate
length_sec = srate_sec*num_samples
length_min = length_sec/60


timespan = np.arange(0, length_sec, srate_sec)

# Plot all signals
plotAll(timespan, meg_data, meg_label, 'Raw Data', 'Time [sec]', 'MEG Data')

# Filter from 0.5-59Hz and Plot
filtered_meg_data = np.zeros_like(meg_data)

for i in range(14):
    filtered_data = amri_sig_filtfft(meg_data[:, i], srate, 0.5, 59)
    filtered_meg_data[:, i] = filtered_data

# Plot the filtered data
plotAll(timespan, filtered_meg_data, meg_label, 'Filtered from 0.5-59Hz', 'Time [sec]', 'MEG Data')


# Isolate the 60-70 sec time frame
start_index = np.argmax(timespan >= 60)
end_index = np.argmax(timespan >= 70)

time_6070 = timespan[start_index:end_index + 1]
filtered_meg_data_6070 = filtered_meg_data[start_index:end_index + 1, :]


# Perform FFT on 60-70 sec window
all_data_6070_FFT = np.zeros((np.shape(filtered_meg_data_6070)[0]/2)+1, np.shape(filtered_meg_data_6070)[1])
for i in range(14):
    data_6070_FFT = np.abs(np.fft.fft(filtered_meg_data_6070[:,i],axis=0))
    all_data_6070_FFT[:,i] = data_6070_FFT

frequency_axis_6070 = np.transpose(np.fft.fftfreq(data_6070_FFT.shape[0], 1/srate))
plotAll(frequency_axis_6070,all_data_6070_FFT,meg_label,'FFT of Data at 60-70 Seconds (No Hanning)', 'Frequency [Hz]', 'Amplitude')


# Apply Hanning window and redo FFT
hanning_window = np.hanning(len(filtered_meg_data_6070))
data_6070_hanning = filtered_meg_data_6070 * hanning_window

all_data_6070_hanning_FFT = np.zeros((np.shape(filtered_meg_data_6070)[0]/2)+1, np.shape(filtered_meg_data_6070)[1])
for i in range(14):
    data_6070_hanning_FFT = np.abs(np.fft.fft(data_6070_hanning,axis=0))
    all_data_6070_hanning_FFT[:,i] = data_6070_hanning_FFT
    
frequency_axis = np.transpose(np.fft.fftfreq(data_6070_hanning_FFT.shape[0], 1/srate))
plotAll(frequency_axis,all_data_6070_hanning_FFT,meg_label,'FFT of Data at 60-70 Seconds (With Hanning)', 'Frequency [Hz]', 'Amplitude')

# Repeat for 60-120 seconds
start_index = np.argmax(timespan >= 60)
end_index = np.argmax(timespan >= 120)

time_60120 = timespan[start_index:end_index + 1]
filtered_meg_data_60120 = filtered_meg_data[start_index:end_index + 1, :]

all_data_60120_FFT = np.zeros((np.shape(filtered_meg_data_60120)[0]/2)+1, np.shape(filtered_meg_data_60120)[1])
for i in range(14):
    data_60120_FFT = np.abs(np.fft.fft(filtered_meg_data_60120,axis=0))
    all_data_60120_FFT[:,i] = data_60120_FFT

frequency_axis_60120 = np.transpose(np.fft.fftfreq(data_60120_FFT.shape[0], 1/srate))
plotAll(frequency_axis_60120,all_data_60120_FFT,meg_label,'FFT of Data at 60-120 Seconds (No Hanning)', 'Frequency [Hz]', 'Amplitude')

hanning_window = np.hanning(len(filtered_meg_data_60120))
data_60120_hanning = filtered_meg_data_60120 * hanning_window

all_data_60120_hanning_FFT = np.zeros((np.shape(filtered_meg_data_60120)[0]/2)+1, np.shape(filtered_meg_data_60120)[1])
for i in range(14):
    data_60120_hanning_FFT = np.abs(np.fft.fft(data_60120_hanning,axis=0))
    all_data_60120_hanning_FFT[:,i] = data_60120_hanning_FFT
    
frequency_axis = np.transpose(np.fft.fftfreq(data_60120_hanning_FFT.shape[0], 1/srate))
plotAll(frequency_axis,all_data_60120_hanning_FFT,meg_label,'FFT of Data at 60-120 Seconds (With Hanning)', 'Frequency [Hz]', 'Amplitude')


# Calculate power spectrum density for 60-70 and 60-120
psd_6070 = np.abs(all_data_6070_FFT)**2 / (len(filtered_meg_data_6070) * srate)
psd_60120 = np.abs(all_data_60120_FFT)**2 / (len(filtered_meg_data_60120) * srate)

plotAll(frequency_axis_6070,psd_6070,meg_label,'PSD of Data at 60-70 Seconds (No Hanning)', 'Frequency [Hz]', 'Power')
plotAll(frequency_axis_60120,psd_60120,meg_label,'PSD of Data at 60-120 Seconds (No Hanning)', 'Frequency [Hz]', 'Power')


# Evaluate PSD between 0-2 minutes
start_index = np.argmax(timespan >= 0)
end_index = np.argmax(timespan >= 120)

time_02 = timespan[start_index:end_index + 1]
filtered_meg_data_02 = filtered_meg_data[start_index:end_index + 1, :]

all_data_02_FFT = np.zeros((np.shape(filtered_meg_data_02)[0]/2)+1, np.shape(filtered_meg_data_02)[1])
for i in range(14):
    data_02_FFT = np.abs(np.fft.fft(filtered_meg_data_02[:,i],axis=0))
    all_data_02_FFT[:,i] = data_02_FFT

frequency_axis_02 = np.transpose(np.fft.fftfreq(data_02_FFT.shape[0], 1/srate))
psd_02 = np.abs(all_data_02_FFT)**2 / (len(filtered_meg_data_02) * srate)
plotAll(frequency_axis_02,psd_02,meg_label,'PSD of Data at 0-2 Minutes', 'Frequency [Hz]', 'Power')



# Evaluate PSD between 4-6 minutes
start_index = np.argmax(timespan >= 240)
end_index = np.argmax(timespan >= 360)

time_46 = timespan[start_index:end_index + 1]
filtered_meg_data_46 = filtered_meg_data[start_index:end_index + 1, :]

all_data_46_FFT = np.zeros((np.shape(filtered_meg_data_46)[0]/2)+1, np.shape(filtered_meg_data_46)[1])
for i in range(14):
    data_46_FFT = np.abs(np.fft.fft(filtered_meg_data_46[:,i],axis=0))
    all_data_46_FFT[:,i] = data_46_FFT

frequency_axis_46 = np.transpose(np.fft.fftfreq(data_46_FFT.shape[0], 1/srate))
psd_46 = np.abs(all_data_46_FFT)**2 / (len(filtered_meg_data_46) * srate)
plotAll(frequency_axis_46,psd_46,meg_label,'PSD of Data at 4-6 Minutes', 'Frequency [Hz]', 'Power')
