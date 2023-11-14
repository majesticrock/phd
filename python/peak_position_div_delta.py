import numpy as np
import matplotlib.pyplot as plt
import lib.continued_fraction as cf
from lib.iterate_containers import naming_scheme
from lib.extract_key import *
import lib.resolvent_peak as rp
# Calculates the resolvent in w^2

prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

Ts = np.array([0.])
Us = np.array([-2.])
Vs = np.array(["50.0", "25.0", "15.0", "10.0", "8.0", "6.0", "4.0", "3.5", "3.0", "2.5", "2.0", "1.5", "1.4", "1.3",
               "1.2", "1.1", "1.0", "0.9", "0.8", "0.7", "0.6", "0.5", 
               "0.45", "0.4", "0.35", "0.3", "0.25", "0.2", "0.15",
               "0.1", "0.07", "0.05", "0.04", "0.03", "0.02", "0.01", 
               "0.007", "0.005", "0.004", "0.003", "0.002", "0.0015", "0.001", 
               "0.0007", "0.0005", "0.0004", "0.0003", "0.0002", "0.00015", "0.0001", 
               "0.00007", "0.00005", "0.00003"])

folder = "data/modes/square/dos_900/"
element_names = ["a", "a+b", "a+ib"]

name_suffix = "phase_SC"

peak_positions = np.zeros(len(Vs))
counter = 0
for name in naming_scheme(Ts, Us, Vs):
    peak_positions[counter] = rp.Peak(f"{folder}{name}", name_suffix).peak_position / extract_key(f"{folder}{name}/resolvent_{name_suffix}.dat.gz", "Total Gap")
    counter += 1

fig, ax = plt.subplots()
cut = -15
v_data = np.log(np.array([float(v) for v in Vs]))
from scipy.optimize import curve_fit
def func(x, a, b):
    return a * x + b

peak_positions = np.log(peak_positions)
popt, pcov = curve_fit(func, v_data[cut:], peak_positions[cut:])
x_lin = np.linspace(np.min(v_data), np.max(v_data), 2000)
ax.plot(x_lin, func(x_lin, *popt), label="Fit")
ax.text(-2, -2, f"$a={popt[0]:.5f}$")
ax.text(-2, -2.7, f"$b={popt[1]:.5f}$")

ax.plot(v_data[cut:], peak_positions[cut:], "X", label="Data Fit")
ax.plot(v_data[:cut], peak_positions[:cut], "o", label="Omitted data")

ax.set_xlabel(r"$\ln(V / t)$")
ax.set_ylabel(r"$\ln(z_0 / \Delta)$")
legend = plt.legend(loc=2)

fig.tight_layout()

import os
plt.savefig(f"python/build/{os.path.basename(__file__).split('.')[0]}.pdf")
plt.show()
