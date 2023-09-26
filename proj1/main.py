import numpy as np
from scipy.io import loadmat
from scipy.fft import rfft, rfftfreq
from scipy.signal import stft
from scipy.integrate import trapezoid
import seaborn as sns
from itertools import cycle 
import matplotlib.pyplot as plt

sns.set()

# Define PlotAll function
def plotAll(x_data, data, dataLabels, figTitle, figXLabel, figYLabel, xlims=None):
    fig = plt.figure(figsize=(12, 14))
    t = fig.add_gridspec(7, 2, wspace=0.05, hspace=0.05)

    colors = plt.cm.plasma(np.linspace(0, 1, 15))
    color_cycle = cycle(colors)

    timePad = 0.0042 * (x_data[-1] - x_data[0])

    if xlims is not None:
        ylow = np.min(data[(xlims[0] <= x_data) & (x_data <= xlims[1])])
        yhigh = np.max(data[(xlims[0] <= x_data) & (x_data <= xlims[1])])
    else:
        ylow = np.min(data)
        yhigh = np.max(data)

    for i in range(14):
        ax = fig.add_subplot(t[i])
        name = dataLabels[i]

        ax.plot(x_data, data[:, i], label=name, color=next(color_cycle))

        ax.set_xlim(x_data[0] - timePad, x_data[-1] + timePad)
        ax.set_ylim(ylow, yhigh)

        if i != 12 and i != 13:
            ax.set_xticklabels([])

        if i % 2 == 1:
            ax.set_yticklabels([])

        leg = ax.legend(handlelength=0, handletextpad=0, fancybox=True, loc='upper right')
        for item in leg.legend_handles:
            item.set_visible(False)

        if xlims is not None:
            ax.set_xlim(xlims)

    fig.text(0.5, 0.04, figXLabel, ha='center')
    fig.text(0.04, 0.5, figYLabel, va='center', rotation='vertical')
    fig.suptitle(figTitle, fontsize=16)
    #plt.show()

def plotSTFT(t, f, Zxx_all, meg_label, figTitle, figXLabel, figYLabel,):
    fig = plt.figure(figsize=(12, 14))
    tile = fig.add_gridspec(7, 2, wspace=0.5, hspace=0.05)
    for i in range(14):
        ax = fig.add_subplot(7, 2, i+1)

        name = meg_label[i]

        mesh = ax.pcolormesh(t, f, Zxx, vmax=vmax, cmap='viridis')
        plt.colorbar(mesh, aspect=10, pad=0.05)
        ax.plot(t, t*0, label=name)
        ax.set_ylim(0, maxFreq)

        if i != 12 and i != 13:
            ax.set_xticklabels([])

        if i % 2 == 1:
            ax.set_yticklabels([])
        
        leg = ax.legend(handlelength=0, handletextpad=0, fancybox=True)
        for item in leg.legend_handles:
            item.set_visible(False)

    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    fig.colorbar(mesh, cax=cbar_ax)

    fig.text(0.5, 0.04, figXLabel, ha='center')
    fig.text(0.04, 0.5, figYLabel, va='center', rotation='vertical')
    fig.suptitle(figTitle, fontsize=16)

## PROBLEM 1 ##
# Import data and find time length in seconds and minutes
output = loadmat('proj1/proj01_data.mat')
meg_data = output['meg_data']
# meg_label from MATLAB array was coming up weird because it was a cell array, so we just rebuilt it
meg_label = ["MLF11", "MZF01", "MRF11", "MLC11", "MZC01", "MRC11", "MLP11", "MZP01", "MRP11", "MLO11", "MZO01", "MRO11", "MLT11", "MRT11"]
meg_locs = output['meg_locs']
srate = output['srate'][0][0]

num_samples = np.shape(meg_data)[0]
srate_sec = 1/srate
length_sec = srate_sec*num_samples
length_min = length_sec/60

timespan = np.arange(0, length_sec, srate_sec)


## PROBLEM 2 ##
# Plot all signals
plotAll(timespan, meg_data, meg_label, 'Raw MEG Data', 'Time [sec]', 'MEG Data [AU]')


## PROBLEM 3 ##
# Filter from 0.5-59Hz and Plot
# Import filtered data from Matlab
filtered_meg_data = np.loadtxt("proj1/filtered_meg_data.csv", delimiter=",", dtype=float)

# Plot the filtered data
plotAll(timespan, filtered_meg_data, meg_label, 'MEG Data, Filtered from 0.5-59Hz', 'Time [sec]', 'MEG Data [AU]')


## PROBLEM 4 ##
# Isolate the 60-70 sec time frame
start_index = np.argmax(timespan >= 60)
end_index = np.argmax(timespan >= 70)

time_6070 = timespan[start_index:end_index + 1]
filtered_meg_data_6070 = filtered_meg_data[start_index:end_index + 1, :]


# Perform FFT on 60-70 sec window
all_data_6070_FFT = np.zeros((1501, 14))
for i in range(14):
    data_6070_FFT = np.abs(rfft(filtered_meg_data_6070[:, i]))
    all_data_6070_FFT[:, i] = data_6070_FFT

frequency_axis_6070 = np.transpose(rfftfreq(filtered_meg_data_6070.shape[0], 1/srate))
plotAll(frequency_axis_6070, all_data_6070_FFT, meg_label, 'FFT of Filtered Data from 60-70 Seconds (No Hanning)', 'Frequency [Hz]', 'Amplitude [AU]', xlims=(0.5, 30))


# Apply Hanning window and redo FFT
hanning_window = np.hanning(len(filtered_meg_data_6070))
all_data_6070_hanning_FFT = np.zeros((1501, 14))
for i in range(14):
    data_6070_hanning = filtered_meg_data_6070[:, i] * hanning_window
    data_6070_hanning_FFT = np.abs(rfft(data_6070_hanning))
    all_data_6070_hanning_FFT[:, i] = data_6070_hanning_FFT
    
plotAll(frequency_axis_6070, all_data_6070_hanning_FFT, meg_label, 'FFT of Filtered Data from 60-70 Seconds (With Hanning)', 'Frequency [Hz]', 'Amplitude [AU]', xlims=(0.5, 30))


## PROBLEM 5 ##
# Repeat for 60-120 seconds
start_index = np.argmax(timespan >= 60)
end_index = np.argmax(timespan >= 120)

time_60120 = timespan[start_index:end_index + 1]
filtered_meg_data_60120 = filtered_meg_data[start_index:end_index + 1, :]

all_data_60120_FFT = np.zeros((9001, 14))
for i in range(14):
    data_60120_FFT = np.abs(rfft(filtered_meg_data_60120[:, i]))
    all_data_60120_FFT[:, i] = data_60120_FFT

frequency_axis_60120 = np.transpose(rfftfreq(filtered_meg_data_60120.shape[0], 1/srate))
plotAll(frequency_axis_60120, all_data_60120_FFT, meg_label, 'FFT of Filtered Data from 60-120 Seconds (No Hanning)', 'Frequency [Hz]', 'Amplitude [AU]', xlims=(0.5, 30))

hanning_window = np.hanning(len(filtered_meg_data_60120))

all_data_60120_hanning_FFT = np.zeros((9001, 14))
for i in range(14):
    data_60120_hanning = filtered_meg_data_60120[:, i] * hanning_window
    data_60120_hanning_FFT = np.abs(rfft(data_60120_hanning, axis=0))
    all_data_60120_hanning_FFT[:, i] = data_60120_hanning_FFT
    

plotAll(frequency_axis_60120, all_data_60120_hanning_FFT, meg_label, 'FFT of Filtered Data from 60-120 Seconds (With Hanning)', 'Frequency [Hz]', 'Amplitude [AU]', xlims=(0.5, 30))


## PROBLEM 6 ##
# Calculate power spectrum density for 60-70 and 60-120
psd_6070 = np.abs(all_data_6070_FFT)**2 / (len(all_data_6070_FFT) * srate)
psd_60120 = np.abs(all_data_60120_FFT)**2 / (len(all_data_60120_FFT) * srate)

plotAll(frequency_axis_6070, psd_6070*(10**6), meg_label, 'Power Spectrum Density of Data from 60-70 Seconds', 'Frequency [Hz]', 'Power [AU^2/Hz*10^-6]', xlims=(0.5, 15))
plotAll(frequency_axis_60120, psd_60120*(10**6), meg_label, 'Power Spectrum Density of Data from 60-120 Seconds', 'Frequency [Hz]', 'Power [AU^2/Hz*10^-6]', xlims=(0.5, 15))


## PROBLEM 7 ##
# Evaluate PSD between 0-2 minutes
start_index = np.argmax(timespan >= 0)
end_index = np.argmax(timespan >= 120)

time_02 = timespan[start_index:end_index + 1]
filtered_meg_data_02 = filtered_meg_data[start_index:end_index + 1, :]

all_data_02_FFT = np.zeros((18001, 14))
for i in range(14):
    data_02_FFT = np.abs(rfft(filtered_meg_data_02[:, i], axis=0))
    all_data_02_FFT[:, i] = data_02_FFT

frequency_axis_02 = np.transpose(rfftfreq(filtered_meg_data_02.shape[0], 1/srate))
psd_02 = np.abs(all_data_02_FFT)**2 / (np.shape(filtered_meg_data_02)[0] * srate)
plotAll(frequency_axis_02, psd_02*(10**7), meg_label, 'Power Spectrum Density of Data from 0-120 Seconds', 'Frequency [Hz]', 'Power [AU^2/Hz*10^-7]', xlims=(0.5, 15))


# Evaluate PSD between 4-6 minutes
start_index = np.argmax(timespan >= 240)
end_index = np.argmax(timespan >= 360)

time_46 = timespan[start_index:end_index + 1]
filtered_meg_data_46 = filtered_meg_data[start_index:end_index + 1, :]

all_data_46_FFT = np.zeros((18001, 14))
for i in range(14):
    data_46_FFT = np.abs(rfft(filtered_meg_data_46[:, i], axis=0))
    all_data_46_FFT[:, i] = data_46_FFT

frequency_axis_46 = np.transpose(rfftfreq(filtered_meg_data_46.shape[0], 1/srate))
psd_46 = np.abs(all_data_46_FFT)**2 / (len(filtered_meg_data_46) * srate)
plotAll(frequency_axis_46, psd_46*(10**6), meg_label, 'Power Spectrum Density of Data from 240-360 Seconds', 'Frequency [Hz]', '[AU^2/Hz*10^-6]', xlims=(0.5, 15))


## PROBLEM 8 ##
maxFreq = 30
length = 30
f, t, Zxx = stft(filtered_meg_data[:length*srate, 0], fs=srate, nperseg=4*srate, noverlap=int(3.5*srate))

Zxx_all = np.zeros((np.size(f), np.size(t), 14))
Zxx_all[:, :, 0] = np.abs(Zxx)

for i in range(1,14):
    f, t, Zxx = stft(filtered_meg_data[:length*srate, i], fs=srate, nperseg=4*srate, noverlap=int(3.5*srate))
    Zxx_all[:, :, i] = np.abs(Zxx)

vmax = np.max(np.abs(Zxx_all)[f<maxFreq, :, :])

plotSTFT(t, f, Zxx_all, meg_label, maxFreq, 'Short Time Fourier Transform from ADD TIME', 'Time [sec])', 'Frequency [Hz]')

## PROBLEM 9 ##
power_all = np.zeros((np.size(t), 14))
for i in range(14):
    power = trapezoid(Zxx_all[(8 <= f) & (f <= 12), :, i], x=f[(8 <= f) & (f <= 12)], axis=0)
    power_all[:,i] = power
    
fig = plt.figure()
plt.plot(t,power_all[:, 10])
plt.xlabel('Time [sec]')
plt.ylabel('Power [AU^2]')
plt.title('8-12Hz Power as Time Series, Channel MZO')


## PROBLEM 10 ##
plotAll(t, power_all, meg_label, '8-12Hz Power as Time Series, All Channels', 'Time [sec]', 'Power [AU^2]')


## PROBLEM 11 ##
average_tf = np.mean(Zxx_all[f<maxFreq, :], axis=2)
vmax = np.max(average_tf)

fig = plt.figure()
plt.pcolormesh(t, f[f<maxFreq], average_tf, vmax=vmax, cmap='viridis')
plt.colorbar()
plt.ylim((0, maxFreq))
plt.xlabel('Time [sec]')
plt.ylabel('Frequency [Hz]')
plt.title('Time-Frequency Analysis,\nAveraged Across All Channels')


## PROBLEM 12 ##
log_average_tf = np.log10(average_tf)
vmax = np.max(log_average_tf)

fig = plt.figure()
plt.pcolormesh(t, f[f<maxFreq], log_average_tf, vmax=vmax, cmap='viridis')
plt.colorbar()
plt.ylim((0, maxFreq))
plt.xlabel('Time [sec]')
plt.ylabel('Frequency [Hz]')
plt.title('Time-Frequency Analysis,\nAveraged Across All Channels (log10 Scale)')

plt.show()
