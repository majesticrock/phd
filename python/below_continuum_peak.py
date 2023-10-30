import numpy as np
import matplotlib.pyplot as plt
import lib.continued_fraction as cf
from lib.iterate_containers import *
from lib.extract_key import * 
# Calculates the resolvent in w^2

prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

realPart = True
both = False

T = 0.0
U = -2.0
V = 0.1
name = f"T={T}/U={U}/V={V}"

folder = "data/modes/square/dos_900/"
fig, ax = plt.subplots()

name_suffix = "phase_SC"

data, data_real, w_lin, res = cf.resolvent_data(f"{folder}{name}", name_suffix, 0.5, 1, number_of_values=50000, imaginary_offset=1e-6, xp_basis=True)
peak_pos = np.argmax(data)
peak_pos_value = w_lin[peak_pos]

from scipy.optimize import fmin
def min_func(x):
    return res.continued_fraction(x + 1e-6j, True).imag
scipy_test = fmin(min_func, peak_pos)

print("Crude:", peak_pos_value)
print("SciPy:", scipy_test[0])
print(min_func(peak_pos_value), min_func(scipy_test[0]))
exit()

data, data_real, w_lin, res = cf.resolvent_data_log_z(f"{folder}{name}", name_suffix, lower_edge=peak_pos_value,
                                                          range=0.5, begin_offset=1e-12,
                                                          number_of_values=2000, imaginary_offset=1e-6, xp_basis=True)
w_lin -= peak_pos_value
plot_data = (data_real)
ax.plot(w_lin, plot_data, label="Data")

from scipy.optimize import curve_fit
def func(x, a, b):
    return a * x + b

popt, pcov = curve_fit(func, w_lin, plot_data)
print(popt[0], " +/- ", np.sqrt(pcov[0][0]))
print(popt[1], " +/- ", np.sqrt(pcov[1][1]))
ax.plot(w_lin, func(w_lin, *popt), "k--", label="Fit")


ax.set_xlabel(r"$z / t$")
ax.set_ylabel(r"Spectral density / a.u.")
#fig.tight_layout()

import os
plt.savefig(f"python/build/{os.path.basename(__file__).split('.')[0]}.pdf")
plt.show()
