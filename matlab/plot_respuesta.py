import numpy as np
import matplotlib.pyplot as plt

from obspy.signal.invsim import paz_to_freq_resp

K=2000
A0=63165

poles = [-0.07405 + 0.07405j, -0.07405 - 0.07405j, -177.72 +177.72j,  -177.72 -177.72j]
zeros = [0.0 + 0.0j, 0.0 + 0.0j]
scale_fac = K*A0

h, f = paz_to_freq_resp(poles, zeros, scale_fac, 0.005, 16384, freq=True)
f = f/(2*np.pi)

plt.figure()
plt.subplot(121)
plt.loglog(f, abs(h))
plt.xlabel('Frequency [Hz]')
plt.ylabel('Amplitude')
plt.axvline(x=1/120, color='r', linestyle='--')
#plt.axvline(x=50, color='r', linestyle='--')
plt.grid(True, which="both", ls="--")


plt.subplot(122)
phase = 2 * np.pi + np.unwrap(np.angle(h))
plt.semilogx(f, phase -2*np.pi)
plt.xlabel('Frequency [Hz]')
plt.ylabel('Phase [radian]')
# ticks and tick labels at multiples of pi
plt.yticks(
    [-np.pi, -np.pi / 2, 0, np.pi / 2, np.pi],
    [r'$-\pi$', r'$-\frac{\pi}{2}$', r'$0$', r'$\frac{\pi}{2}$', r'$\pi$'])
plt.ylim(-np.pi-0.2,np.pi + 0.2)
plt.grid(True, which="both", ls="--")
# title, centered above both subplots
plt.suptitle('Frequency Response of 151B Seismometer')
# make more room in between subplots for the ylabel of right plot
plt.subplots_adjust(wspace=0.3)
plt.show()

