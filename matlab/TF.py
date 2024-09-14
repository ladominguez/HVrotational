import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
K = 2000
A0 = 63165
def H(s):
    N = K * A0 * (s**2)
    D = (s**2 + 0.14807*s + 0.010966) * (s**2 + 355.38*s + 63165)
    return N / D
freq = np.logspace(-2, 2, 100) 
s = 1j * freq * 2 * np.pi
H_s = H(s)
M = np.abs(H_s)
ph = np.angle(H_s, deg = True)
dM = np.gradient(20 * np.log10(M), freq)
umbral = np.std(dM) * 2 
change_indices = np.where(np.abs(np.gradient(dM)) > umbral)[0]
lineas_rojas_x = freq[change_indices]
plt.figure(figsize = (12, 6))
# Magnitud.
plt.subplot(2, 1, 1)
plt.semilogx(freq, 20 * np.log10(M), color = 'green')
for x in lineas_rojas_x:
    plt.axvline(x = x, color = 'red', linestyle = '--')
    plt.text(x, 20 * np.log10(np.min(M)), f'{x:.4f}', color = 'red', rotation = 90, verticalalignment = 'bottom')
plt.title('Función de Transferencia')
plt.ylabel('Magnitud (dB)')
plt.grid(which = 'both', linestyle = '--')
# Fase.
plt.subplot(2, 1, 2)
plt.semilogx(freq, ph, color = 'blue')
for x in lineas_rojas_x:
    plt.axvline(x = x, color = 'red', linestyle = '--')
    plt.text(x, np.min(ph), f'{x:.4f}', color = 'red', rotation = 90, verticalalignment = 'bottom')
plt.ylabel('Fase (grados)')
plt.xlabel('Frecuencia (Hz)')
plt.grid(which = 'both', linestyle = '--')
plt.show()
