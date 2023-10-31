import numpy

## 1. How long is the recording in seconds? 

## 2. For each channel, evaluate its mean and standard deviation? 

## 3. Evaluate the correlation between every pair of different channels. 

## 4. For each channel, detect the neuronal spikes. Hint: this is an open-ended question. 
## You need to observe the data first and come up with your own criteria for a spike. 
## Typically, a spike will show a peak that deviates from the baseline (often just the mean) by a 
## notable extent (e.g., multiple times of a standard deviation). The peak may be either positive or negative. Exercise your own judgment. 

## 5. For each spike you have detected in 4, choose a short window around the peak of the spike. 
## The window is often very short. You may explore the data and exercise your own judgment. 
## For example, a window may start from 0.5 ms before the peak and 0.5 ms after the peak. 
## The signal within such a window defines the shape of each spike you have detected. 

## 6. Find a way to cluster the spike shape using k-means clustering. Each cluster corresponds to one neuron. 
## Hint: you may define features of each spike using PCA or SVD, or other ways you would come up with. This is again open-ended. 

## 7. Justify the number of clusters (or neurons) Calculate the centroid of each cluster. 
## Use this centroid to obtain the average spike shape of each neuron. 

## 8. Repeat the analysis for each channel separately. 

## 9. Calculate the firing rate of each neuron detected. That is to count how many times each neuron fire within a 10-s moving window. 

## 10. Evaluate the temporal fluctuation of each neuron's firing rate. 

## 11. Explore the relationship between the firing rate and the stimulation delivered to the peripheral nerve. 