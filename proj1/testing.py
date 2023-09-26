import numpy as np
from scipy.fft import rfft, irfft, rfftfreq
import pylab as plt

time   = np.arange(0,10,1/300)
signal = np.cos(5*np.pi*time) + np.cos(7*np.pi*time)

freq = rfftfreq(np.size(signal), 1/300)
fftsig = rfft(signal)


plt.subplot(211)
plt.plot(time,signal)
plt.subplot(212)
plt.plot(freq,fftsig)
plt.xlim((0, 10))
plt.show()
