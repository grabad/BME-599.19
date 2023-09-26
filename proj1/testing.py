import numpy as np
from scipy.io import loadmat
from scipy.fft import rfft, rfftfreq
from scipy.signal import stft
import seaborn as sns
from itertools import cycle 
import matplotlib.pyplot as plt

sns.set()

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

maxFreq = 30

fig = plt.figure(figsize=(12, 14))
t = fig.add_gridspec(7, 2, wspace=0.5, hspace=0.05)
for i in range(14):
    ax = fig.add_subplot(7, 2, i+1)

    name = meg_label[i]

    f, t, Zxx = stft(filtered_meg_data[:3000, i], fs=srate, nperseg=4*srate, noverlap=int(3.5*srate))
    vmax = np.max(np.abs(Zxx)[f<maxFreq,:])
    
    ax.pcolormesh(t, f, np.abs(Zxx), vmax=vmax)
    ax.plot(t, t*0, label=name)
    ax.set_ylim(0, maxFreq)

    if i != 12 and i != 13:
        ax.set_xticklabels([])

    if i % 2 == 1:
        ax.set_yticklabels([])
    
    leg = ax.legend(handlelength=0, handletextpad=0, fancybox=True)
    for item in leg.legend_handles:
        item.set_visible(False)

    # if i == 5:
    #     fig.colorbar()

fig.text(0.5, 0.04, 'Time (s)', ha='center')
fig.text(0.04, 0.5, 'Frequency (Hz)', va='center', rotation='vertical')
fig.suptitle('STFT', fontsize=16)
plt.show()