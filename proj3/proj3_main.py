import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat
import scipy.signal as ss
import seaborn as sns
from sklearn.cluster import KMeans
from itertools import chain

sns.set()

## 1. How long is the recording in seconds? 
raw_data = loadmat('proj3/data.mat')

data = raw_data['data']
fs = raw_data['fs'][0][0]
stim = raw_data['stim']
data_length = np.size(data, axis=1)
time_length = data_length/fs

print(f'Data is {time_length} seconds long\n')

## 2. For each channel, evaluate its mean and standard deviation? 

means = np.mean(data, axis=1)
stdevs = np.std(data, axis=1)

for i in range(5):
    print(f'Channel {i+1}: Mean={means[i]}, STD={stdevs[i]}')
print()

## 3. Evaluate the correlation between every pair of different channels. 
pairs =  np.array([[0, 1], 
          [0, 2], 
          [0, 3], 
          [0, 4], 
          [1, 2], 
          [1, 3], 
          [1, 4], 
          [2, 3], 
          [2, 4], 
          [3, 4]])

coefs = np.zeros(10)

for i in range(10):
    coefs[i] = np.corrcoef(data[pairs[i, 0]], data[pairs[i, 1]])[0, 1]
    print(f'Correlation between Channel {pairs[i, 0]+1} and {pairs[i, 1]+1}: {coefs[i]}')
print()

## 4. For each channel, detect the neuronal spikes. Hint: this is an open-ended question. 
## You need to observe the data first and come up with your own criteria for a spike. 
## Typically, a spike will show a peak that deviates from the baseline (often just the mean) by a 
## notable extent (e.g., multiple times of a standard deviation). The peak may be either positive or negative. Exercise your own judgment. 
def find_spikes(data, mean, stdev, fs, window):
    peaks_p, _ = ss.find_peaks(data, height=(mean+2*stdev), distance=int(fs*window))
    peaks_n, _ = ss.find_peaks(-data, height=(-mean+2*stdev), distance=int(fs*window))
    peaks = np.concatenate((peaks_p, peaks_n))
    peaks.sort()

    return peaks

window = 0.005
peaks = [find_spikes(data[0, :], means[0], stdevs[0], fs, window), 
        find_spikes(data[1, :], means[1], stdevs[1], fs, window), 
        find_spikes(data[2, :], means[2], stdevs[2], fs, window), 
        find_spikes(data[3, :], means[3], stdevs[3], fs, window), 
        find_spikes(data[4, :], means[4], stdevs[4], fs, window)]

## 5. For each spike you have detected in 4, choose a short window around the peak of the spike. 
## The window is often very short. You may explore the data and exercise your own judgment. 
## For example, a window may start from 0.5 ms before the peak and 0.5 ms after the peak. 
## The signal within such a window defines the shape of each spike you have detected. 
def spike_shape(data, spike_index, fs, window):
    start = int(spike_index - fs*window) - 1
    if start < 0: start = 0

    stop = int(spike_index + fs*window)
    if stop > len(data): stop = len(data)

    spike = data[start:stop]

    return spike

spike_0 = spike_shape(data[0, :], peaks[0][0], fs, window)

## 6. Find a way to cluster the spike shape using k-means clustering. Each cluster corresponds to one neuron. 
## Hint: you may define features of each spike using PCA or SVD, or other ways you would come up with. This is again open-ended. 
spike_means = np.zeros(len(list(chain(*peaks))))
spike_vars = np.zeros(len(list(chain(*peaks))))
count = 0

for row_index in range(5):
    for peak_index in peaks[row_index]:
        spike = spike_shape(data[row_index, :], peak_index, fs, window)

        spike_means[count] = np.mean(spike)
        spike_vars[count] = np.var(spike)

        count = count + 1

cluster_data = zip(spike_means, spike_vars)

## 7. Justify the number of clusters (or neurons) Calculate the centroid of each cluster. 
## Use this centroid to obtain the average spike shape of each neuron. 

## 8. Repeat the analysis for each channel separately. 

## 9. Calculate the firing rate of each neuron detected. That is to count how many times each neuron fire within a 10-s moving window. 

## 10. Evaluate the temporal fluctuation of each neuron's firing rate. 

## 11. Explore the relationship between the firing rate and the stimulation delivered to the peripheral nerve. 
