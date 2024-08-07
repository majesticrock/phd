import numpy as np
import matplotlib.pyplot as plt
import __path_appender as __ap
__ap.append()
import continued_fraction as cf

T = 0.
U = -2.5
V = 0.0
name = f"T={T}/U={U}/V={V}"

use_xp = True
folder = "data/modes/square/dos_3k_SC/"
name_suffix = "higgs_CDW"
fig, ax = plt.subplots()

plot_lower_lim = 0.08
plot_upper_lim = plot_lower_lim + 0.2
data_imag, data, w_lin, res = cf.resolvent_data(f"{folder}{name}", name_suffix, plot_lower_lim, plot_upper_lim, xp_basis=use_xp, imaginary_offset=0)

from scipy.optimize import curve_fit
def func_ln(x, a, b):
    return a * x + b

try:
    w_log = np.log((w_lin.real))
    ax.plot(w_log, np.log(data), "-", label="Data")
    popt, pcov = curve_fit(func_ln, w_log, np.log(data))
    ax.plot(w_log, func_ln(w_log, *popt), "--", label=r"$a \ln( z ) + b$")
    ax.set_xlabel(r"$\ln(z)$")
    ax.set_ylabel(r"$\ln(\Re G^\mathrm{ret}(z))$")
    ax.text(0.05, 0.35, f"$a={popt[0]:.5f}$", transform = ax.transAxes)
    ax.text(0.05, 0.3, f"$b={popt[1]:.5f}$", transform = ax.transAxes)
except RuntimeError:
    print("Could not estimate curve_fit")
except ValueError:
    print("Value")

ax.legend()
fig.tight_layout()

import os
plt.savefig(f"python/build/{os.path.basename(__file__).split('.')[0]}_U={U}.pdf")
plt.show()