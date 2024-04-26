import numpy as np
import matplotlib.pyplot as pltwe
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib.continued_fraction as cf
from lib.iterate_containers import naming_scheme
from lib.extract_key import *
import scipy.optimize as opt
import lib.resolvent_peak as rp
# Calculates the resolvent in w^2

prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

T = 0.
Us = np.array([3.7, 3.725, 3.75, 3.8, 3.825, 3.85, 3.875, 3.9, 3.91,  # 9
                3.925, 3.94, 3.95, 3.96, 3.97, 3.975, 3.985, 3.99, 3.995, # 9
                3.997, 3.9985, 3.999, 3.9995, # 4
                #4.0, # 1
                4.0005, 4.001, 4.0015, 4.003, # 4
                4.005, 4.01, 4.015, 4.025, 4.03, 4.04, 4.05, 4.06, 4.075, # 9
                4.09, 4.1, 4.125, 4.15, 4.175, 4.2, 4.25, 4.275, 4.3 # 9
    ])
V = 1.

folder = "data/modes/square/dos_3k/"
colors = ["orange", "purple", "green",  "black"]

name_suffix = "CDW"
fig, ax = plt.subplots()

name = f"T={T}/U=4.0/V={V}"
zero_peak = rp.Peak(f"{folder}{name}", name_suffix, (2.45, 2.5))
zero_value = zero_peak.improved_peak_position()["x"]

peak_positions = np.zeros(len(Us))
counter = 0
for name in naming_scheme(T, Us, V):
    peak = rp.Peak(f"{folder}{name}", name_suffix, (2.4, 2.85))
    peak_positions[counter] = np.log(peak.improved_peak_position()["x"] - zero_value)
    counter += 1

cut = -int(len(Us) / 2)
u_data = np.log(np.abs(np.array([float(u) for u in Us]) - 4 * V))

from scipy.optimize import curve_fit
def func(x, a, b):
    return a * x + b

popt, pcov = curve_fit(func, u_data[cut:], peak_positions[cut:])
ax.text(0.05, 0.35, f"$a_\\mathrm{{AFM}}={popt[0]:.5f}$", transform = ax.transAxes)
ax.text(0.05, 0.3, f"$b_\\mathrm{{AFM}}={popt[1]:.5f}$", transform = ax.transAxes)

x_lin = np.linspace(np.min(u_data), np.max(u_data), 2000)
ax.plot(x_lin, func(x_lin, *popt), color=colors[0])
    
cut2 = int(len(Us) / 2)
popt, pcov = curve_fit(func, u_data[:cut2], peak_positions[:cut2])
ax.text(0.05, 0.5, f"$a_\\mathrm{{CDW}}={popt[0]:.5f}$", transform = ax.transAxes)
ax.text(0.05, 0.45, f"$b_\\mathrm{{CDW}}={popt[1]:.5f}$", transform = ax.transAxes)

ax.plot(x_lin, func(x_lin, *popt), color=colors[1]) 
ax.plot(u_data[cut:], peak_positions[cut:], "X", label=f"{name_suffix} peak in AFM", color=colors[0])
ax.plot(u_data[:cut2], peak_positions[:cut2], "X", label=f"{name_suffix} peak in CDW", color=colors[1])

ax.set_xlabel(r"$\ln((U - 4V) / t)$")
ax.set_ylabel(r"$\ln((z_0 - z_0(U=4V)) / t)$")
legend = plt.legend()

fig.tight_layout()

import os
plt.savefig(f"python/build/{os.path.basename(__file__).split('.')[0]}.pdf")
plt.show()
