import numpy as np
import matplotlib.pyplot as plt
from lib.ez_fit import *
import lib.resolvent_peak as rp

T = 0.0
U = -2.5
V = -0.1
name = f"T={T}/U={U}/V={V}"

N = np.array([400, 500, 600, 750, 900, 1200, 1500, 2000, 3000, 6000])
peak_positions = np.array([0.0259, 0.02286, 0.02114, 0.01888, 0.017246, 0.015029, 0.01357, 0.01183, 0.00969, 0.00728])

name_suffix = "phase_SC"
peak_fit_params = {"range": 1e-4, "begin_offset": 1e-10, "imaginary_offset": 1e-7, "peak_position_tol": 1e-14}

for i, n in enumerate(N):
    folder = f"data/modes/square/dos_{n}/"
    peak_positions[i], weight = rp.find_weight(f"{folder}{name}", name_suffix, **peak_fit_params)

fig, ax = plt.subplots()

ax.plot(1/N, peak_positions, "o", label="Data")
def func(x, a, b):
    return a*np.sqrt(x) + b
popt, pcov = ez_general_fit(1/N, peak_positions, func, label="$a \\sqrt{x} + b$")

ax.text(0.00175, 0.0175, f"$a = {popt[0]:.5f}$")
ax.text(0.00175, 0.0155, f"$b = {popt[1]:.5f}$")

ax.set_xlabel("$1/N_\\gamma$")
ax.set_ylabel("$\\omega_0$")

ax.legend()
plt.tight_layout()
plt.show()