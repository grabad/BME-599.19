import numpy as np
from scipy.io import loadmat
from scipy.fft import rfft, rfftfreq
from scipy.signal import stft
from scipy.integrate import trapezoid
import seaborn as sns
from itertools import cycle 
import matplotlib.pyplot as plt

sns.set()

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

def plotSTFT(t, f, Zxx_all, maxFreq, meg_label, figTitle, figXLabel, figYLabel, ):
    fig = plt.figure(figsize=(12, 14))
    tile = fig.add_gridspec(7, 2, wspace=0.5, hspace=0.05)
    for i in range(14):
        ax = fig.add_subplot(7, 2, i+1)
        name = meg_label[i]

        Zxx = Zxx_all[:, :, i]
        vmax = np.max(Zxx[f<maxFreq, :])

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

    fig.text(0.5, 0.04, figXLabel, ha='center')
    fig.text(0.04, 0.5, figYLabel, va='center', rotation='vertical')
    fig.suptitle(figTitle, fontsize=16)

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
filtered_meg_data = np.loadtxt("proj1/filtered_meg_data.csv", delimiter=",", dtype=float)

# plotAll(timespan, filtered_meg_data, meg_label, 'Filtered from 0.5-59Hz', 'Time (s)', 'MEG Data')

maxFreq = 60
length = 2400
f, t, Zxx = stft(filtered_meg_data[:length*srate, 0], fs=srate, nperseg=4*srate, noverlap=int(3.5*srate))

Zxx_all = np.zeros((np.size(f), np.size(t), 14))
Zxx_all[:, :, 0] = np.abs(Zxx)

for i in range(1,14):
    f, t, Zxx = stft(filtered_meg_data[:length*srate, i], fs=srate, nperseg=4*srate, noverlap=int(3.5*srate))
    Zxx_all[:, :, i] = np.abs(Zxx)

vmax = np.max(np.abs(Zxx_all)[f<maxFreq, :, :])

plotSTFT(t, f, Zxx_all, maxFreq, meg_label, 'Short Time Fourier Transform, All Channels', 'Time [sec]', 'Frequency [Hz]')

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