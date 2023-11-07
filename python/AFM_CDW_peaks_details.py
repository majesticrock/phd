import numpy as np
import matplotlib.pyplot as plt
import lib.continued_fraction as cf
from lib.iterate_containers import iterate_containers
from lib.extract_key import *
import scipy.optimize as opt
# Calculates the resolvent in w^2

prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

T = 0.
Us = np.array([3.7, 3.75, 3.8, 3.85, 3.9, 3.925, 3.95, 3.975, 3.99, 3.995, # 10 values
                    3.9975, 3.9985, 3.999, 3.9995, 3.99975, 3.99985, 3.9999, # 7
                    4.0001, 4.00015, 4.00025, 4.0005, 4.001, 4.0015, 4.0025, # 7
                    4.005, 4.01, 4.025, 4.05, 4.075, 4.1, 4.15, 4.2, 4.25, 4.3]) # 10
V = 1.

folder = "data/modes/square/dos_900/"
colors = ["orange", "purple", "green",  "black"]

name_suffix = "CDW"
fig, ax = plt.subplots()

name = f"T={T}/U=4.0/V={V}"
data, data_real, w_lin, res = cf.resolvent_data(f"{folder}{name}", name_suffix, lower_edge=2.45, upper_edge=2.5, number_of_values=20000, imaginary_offset=1e-6, xp_basis=True, messages=False)
zero_value = w_lin[np.argmax(data)]

peak_positions = np.zeros(len(Us))
counter = 0
for U in Us:
    name = f"T={T}/U={U}/V={V}"
    data, data_real, w_lin, res = cf.resolvent_data(f"{folder}{name}", name_suffix, lower_edge=2.4, upper_edge=2.85, number_of_values=20000, imaginary_offset=1e-6, xp_basis=True, messages=False)

    def min_func(x):
        return res.continued_fraction(x + 1e-6j, True).imag
    scipy_result = opt.fmin_l_bfgs_b(min_func, w_lin[np.argmax(data)], bounds=[(2.4, 2.85)], approx_grad=True, epsilon=1e-10)
    #print(scipy_result[0][0], w_lin[np.argmax(data)], " -> ", scipy_result[0][0] - w_lin[np.argmax(data)])
    gap_value = 2 * extract_key(f"{folder}{name}/resolvent_higgs_{name_suffix}.dat.gz", "Total Gap")
    peak_positions[counter] = (scipy_result[0][0] - zero_value)
    counter += 1

cut = -int(len(Us) / 2)
u_data = np.log(np.abs(np.array([float(u) for u in Us]) - 4 * V))

from scipy.optimize import curve_fit
def func(x, a, b):
    return a * x + b

peak_positions = np.log(peak_positions)
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
