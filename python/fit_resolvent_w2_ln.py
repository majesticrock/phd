import numpy as np
import matplotlib.pyplot as plt
import lib.continued_fraction as cf

T = 0.
U = -2.0
V = -0.1
name = f"T={T}/U={U}_V={V}"

use_xp = False
if use_xp:
    folder = "data/modes/square/L=70/"
else:
    folder = "data/modes/square/dos_900/"
name_suffix = "phase_SC"
fig, ax = plt.subplots()

data_imag, data, w_lin, res = cf.resolvent_data(f"{folder}{name}", name_suffix, 0, 1, xp_basis=use_xp)
plot_lower_lim = w_lin[np.argmax(data)] + 0.05
plot_upper_lim = plot_lower_lim + 0.2
data_imag, data, w_lin, res = cf.resolvent_data(f"{folder}{name}", name_suffix, plot_lower_lim, plot_upper_lim, xp_basis=use_xp)

from scipy.optimize import curve_fit
def func_ln(x, a, b):
    return a * x + b

try:
    w_log = np.log(w_lin.real)
    ax.plot(w_log, np.log(data), "-", label="Real")
    popt, pcov = curve_fit(func_ln, w_log, np.log(data))
    ax.plot(w_log, func_ln(w_log, *popt), "--", label=r"$a \ln( z ) + b$")
    ax.set_xlabel(r"$\ln(z)$")
    ax.set_ylabel(r"$\ln(g(z))$")
    ax.text(0.05, 0.35, f"$a={popt[0]}$", transform = ax.transAxes)
    ax.text(0.05, 0.3, f"$b={popt[1]}$", transform = ax.transAxes)
except RuntimeError:
    print("Could not estimate curve_fit")
except ValueError:
    print("Value")

ax.legend()
fig.tight_layout()

import os
plt.savefig(f"python/build/{os.path.basename(__file__).split('.')[0]}_U={U}.pdf")
plt.show()