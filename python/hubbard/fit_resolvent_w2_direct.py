import numpy as np
import matplotlib.pyplot as plt
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib.continued_fraction as cf

T = 0.
U = -3.0
V = -0.1
name = f"T={T}/U={U}/V={V}"

use_xp = False
if use_xp:
    folder = "data/modes/square/L=70/"
else:
    folder = "data/modes/square/momentum_L=20/"
name_suffix = "phase_SC"
fig, ax = plt.subplots()

data_imag, data, w_lin = cf.resolvent_data(f"{folder}{name}", name_suffix, 0, 1, xp_basis=use_xp)
plot_lower_lim = w_lin[np.argmax(data)] + 0.25
plot_upper_lim = plot_lower_lim + 0.2
data_imag, data, w_lin = cf.resolvent_data(f"{folder}{name}", name_suffix, plot_lower_lim, plot_upper_lim, xp_basis=use_xp)

from scipy.optimize import curve_fit
def func(x, a, b):
    return np.sign(x) * a * np.abs(x)**b

try:
    ax.plot(w_lin, data, "-", label="Real")
    popt, pcov = curve_fit(func, w_lin, data)
    print(f"{popt[0]} +/- {pcov[0][0]}")
    print(popt)
    ax.plot(w_lin, func(w_lin, *popt), "--", label=r"$a(x - x_0)^b$")
    ax.set_xlabel(r"$z$")
    ax.set_ylabel(r"$\Re [G]$")
    ax.text(0.05, 0.35, f"$a={popt[0]}$", transform = ax.transAxes)
    ax.text(0.05, 0.3, f"$b={popt[1]}$", transform = ax.transAxes)
    ax.set_yscale("log")
    ax.set_xscale("log")
except RuntimeError:
    print("Could not estimate curve_fit")
except ValueError:
    print("Value")

ax.legend()
fig.tight_layout()

import os
plt.savefig(f"python/build/{os.path.basename(__file__).split('.')[0]}_U={U}.pdf")
plt.show()