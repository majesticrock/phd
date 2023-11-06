import numpy as np
import matplotlib.pyplot as plt
import lib.continued_fraction as cf
from lib.iterate_containers import iterate_containers
from lib.extract_key import *
# Calculates the resolvent in w^2

prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

Ts = np.array([0.])
Us = np.array([3.75, 3.8, 3.85, 3.9, 3.925, 3.95, 3.975, 3.99, 3.995, 3.999,
               4.001, 4.005, 4.01, 4.025, 4.05, 4.075, 4.1, 4.15, 4.2, 4.25])
Vs = np.array([1.])

folder = "data/modes/square/dos_900/"
element_names = ["a", "a+b", "a+ib"]
colors = ["orange", "purple", "green",  "black"]

name_suffix = "CDW"
fig, ax = plt.subplots()

peak_positions = np.zeros(len(Us))
counter = 0
for T, U, V in iterate_containers(Ts, Us, Vs):
    name = f"T={T}/U={U}/V={V}"
    data, data_real, w_lin, res = cf.resolvent_data(f"{folder}{name}", name_suffix, lower_edge=2, upper_edge=3, number_of_values=20000, imaginary_offset=1e-6, xp_basis=True, messages=False)

    gap_value = 2 * extract_key(f"{folder}{name}/resolvent_higgs_{name_suffix}.dat.gz", "Total Gap")
    peak_positions[counter] = w_lin[np.argmax(data)] / gap_value
    counter += 1

cut = -int(len(Us) / 2)
u_data = (np.abs(np.array([float(u) for u in Us]) - 4 * Vs[0]))
ax.set_xscale("log")
from scipy.optimize import curve_fit
def func(x, a, b):
    return a * x + b

peak_positions = np.log(peak_positions)
popt, pcov = curve_fit(func, u_data[cut:], peak_positions[cut:])
x_lin = np.linspace(np.min(u_data), np.max(u_data), 200)
ax.plot(x_lin, func(x_lin, *popt), color=colors[0])
print(name_suffix, popt)
    
cut2 = int(len(Us) / 2)
popt, pcov = curve_fit(func, u_data[:cut2], peak_positions[:cut2])
ax.plot(x_lin, func(x_lin, *popt), color=colors[1])
print(name_suffix, popt)
    
ax.plot(u_data[cut:], peak_positions[cut:], "X", label=f"{name_suffix} peak in AFM", color=colors[0])
ax.plot(u_data[:cut2], peak_positions[:cut2], "X", label=f"{name_suffix} peak in CDW", color=colors[1])

ax.set_xlabel(r"$\ln(U - 4V) / t$")
ax.set_ylabel(r"$z_0 / t$")
legend = plt.legend()

fig.tight_layout()


import os
plt.savefig(f"python/build/{os.path.basename(__file__).split('.')[0]}.pdf")
plt.show()
