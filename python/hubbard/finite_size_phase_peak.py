import numpy as np
import matplotlib.pyplot as plt
import __path_appender as __ap
__ap.append()
from ez_fit import *
import resolvent_peak as rp

T = 0.0
U = -2.5
V = -0.1
name = f"T={T}/U={U}/V={V}"

N = np.array([400, 500, 600, 750, 900, 1200, 1500, 2000, 3000, 6000])
peak_positions = np.array([0.0259, 0.02286, 0.02114, 0.01888, 0.017246, 0.015029, 0.01357, 0.01183, 0.00969, 0.00728])

name_suffix = "phase_SC"
peak_fit_params = {"range": 1e-4, "begin_offset": 1e-10, "imaginary_offset": 1e-7, "peak_position_tol": 1e-14}
lattice_folders = ["square", "cube"]

fig, ax = plt.subplots()

ax.set_ylim(0, 0.03)

for j, lattice_folder in enumerate(lattice_folders):
    for i, n in enumerate(N):
        folder = f"data/modes/{lattice_folder}/dos_{n}/"
        peak = rp.Peak(f"{folder}{name}", name_suffix, imaginary_offset=1e-7)
        peak_pos_value = np.copy(peak.peak_position)
        peak_result = peak.improved_peak_position(xtol=1e-14)
        # only an issue if the difference is too large;
        if not peak_result["success"]:
            print("We might not have found the peak for data_folder!\nWe found ", peak_pos_value, " and\n", peak_result)
        
        peak_positions[i] = peak_result["x"]

    ax.plot(1/N, peak_positions, "o", label=f"Data {lattice_folder}")
    if(j==0):
        def func(x, a, b):
            return a*np.sqrt(x) + b
        popt, pcov = ez_general_fit(1/N, peak_positions, func, label="$a \\sqrt{x} + b$")
        ax.text(0.00175, 0.0175, f"$a = {popt[0]:.5f}$")
        ax.text(0.00175, 0.0155 , f"$b = {popt[1]:.5f}$")

ax.set_xlabel("$1/N_\\gamma$")
ax.set_ylabel("$\\omega_0$")

ax.legend()
plt.tight_layout()
plt.show()