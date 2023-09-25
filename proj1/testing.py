import numpy as np
from scipy.fft import rfft, irfft, rfftfreq

def fftFilt(data, lowCO=None, highCO=None):
    freqData = rfft(data)
    freqRange = rfftfreq(data.size, 1/3000)

    if lowCO is not None:
        freqData[abs(freqRange) < lowCO] = 0
        
    if highCO is not None:
        freqData[abs(freqRange) > highCO] = 0
    
    filteredData = irfft(freqData)
    
    return (filteredData, freqRange)

time   = np.linspace(0,10,1/300)
signal = np.cos(5*np.pi*time) + np.cos(7*np.pi*time)

cut_signal, freqRange = fftFilt(signal, lowCO=6)

import pylab as plt
plt.subplot(211)
plt.plot(time,signal)
plt.subplot(212)
plt.plot(time,cut_signal)
plt.show()
