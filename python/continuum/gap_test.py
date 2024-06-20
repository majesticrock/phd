import matplotlib.pyplot as plt
import numpy as np
import gzip

MEV_FACTOR = 1e3
plotPhase = False

fig, ax = plt.subplots()

with gzip.open(f"data/continuum/test/gap.dat.gz", 'rt') as f_open:
    M = np.loadtxt(f_open)
M[0] = M[0] - M[0][int(0.5 * len(M[0])) - 1]

if plotPhase:
    delta_abs = np.sqrt(M[1]**2 + M[2]**2)
    delta_phase = np.arctan2(M[2], M[1])
    M[1] = delta_abs
    ax.plot(M[0], MEV_FACTOR * delta_abs, "-", label=r"Abs SC")
    ax2 = ax.twinx()
    ax2.plot(M[0], delta_phase / np.pi, "--", color="C3", label=r"Phase SC")
    ax2.set_ylabel(r"arg[$\Delta$] / $\pi$")
else:
    ax.plot(M[0], MEV_FACTOR * M[1], "-", label=r"$\Delta_\mathrm{SC}$")
ax.plot(M[0], MEV_FACTOR * M[3], "-", label=r"$\Delta_n$")

from scipy.optimize import curve_fit
use_points = .6
fit_data_x = M[0][int(use_points * len(M[0])):]
fit_data_y = MEV_FACTOR * M[1][int(use_points * len(M[1])):]

def func(x, a, b):
    return a / x**b
try:
    popt, pcov = curve_fit(func, fit_data_x, fit_data_y)
    print(popt)
    ax.plot(fit_data_x, func(fit_data_x, *popt), ls="--", color="C3", label="Fit")
except RuntimeError:
    print("Fit was not possible")

ax.set_xlabel(r"$k - k_\mathrm{F} [\sqrt{\mathrm{eV}}]$")
ax.set_ylabel(r"$\Delta [\mathrm{meV}]$")
fig.legend()
fig.tight_layout()

plt.show()