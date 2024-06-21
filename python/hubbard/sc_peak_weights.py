import numpy as np
import matplotlib.pyplot as plt
import __path_appender as __ap
__ap.append()
from iterate_containers import *
from extract_key import *
import resolvent_peak as rp

prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

Ts = np.array([0.])
Us = np.array([-2.])
Vs = np.array(["50.0", "25.0", "15.0", "10.0", "8.0", "6.0", "4.0", "3.5", "3.0", "2.5", "2.0", "1.5", 
                "1.4", "1.3", "1.2", "1.1", "1.0", "0.9", "0.8", "0.7", "0.6", "0.5", 
                "0.45", "0.4", "0.35", "0.3", "0.25", "0.2", "0.15", "0.1", "0.07", 
                "0.05", "0.04", "0.03", "0.02", "0.01",
                "0.007", "0.005", "0.004", "0.003", "0.002", "0.0015", "0.001", 
                "0.0007", "0.0005", "0.0004", "0.0003", "0.0002", "0.00015", "0.0001", 
                "0.00007", "0.00005", "0.00003"
                ])

folder = "data/modes/square/dos_900/"
name_suffix = "phase_SC"

weights = np.zeros(len(Vs))
counter = 0

#fig_check, ax_check = plt.subplots()
for T, U, V in iterate_containers(Ts, Us, Vs):
    name = f"T={T}/U={U}/V={V}"
    lower = float(V)
    upper = 8 * float(V) + 2
    
    peak = rp.Peak(f"{folder}{name}", name_suffix, (lower, upper))
    peak_pos_value = np.copy(peak.peak_position)
    peak_result = peak.improved_peak_position(xtol=1e-13)
    # only an issue if the difference is too large;
    if not peak_result["success"]:
        print(f"We might not have found the peak for V={V}!\nWe found ", peak_pos_value, " and\n", peak_result)
    peak_pos_value = peak_result["x"]
    popt, pcov, w_log, y_data = peak.fit_real_part(0.01, 1e-8)
    weights[counter] = popt[1]
    counter += 1
    
    #line, = ax_check.plot(w_log, y_data, ls="-", label=f"V={V}")
    #ax_check.plot(w_log, rp.linear_function(w_log, *popt), ls="--", linewidth=4, color=line.get_color())

v_data = np.log(np.array([float(v) for v in Vs]))
cut = -18
fig, ax = plt.subplots()

ax.plot(v_data[cut:], weights[cut:], "X", label="Fitted data")
ax.plot(v_data[:cut], weights[:cut], "o", label="Omitted data")

from scipy.optimize import curve_fit

popt, pcov = curve_fit(rp.linear_function, v_data[cut:], weights[cut:])
v_lin = np.linspace(v_data.min(), v_data.max(), 500)
ax.plot(v_lin, rp.linear_function(v_lin, *popt), label="Fit")

ax.text(0.05, 0.35, f"$a={popt[0]:.5f}$", transform = ax.transAxes)
ax.text(0.05, 0.3, f"$b={popt[1]:.5f}$", transform = ax.transAxes)

ax.set_xlabel(r"$\ln(V / t)$")
ax.set_ylabel(r"$\ln(w_0 \cdot t)$")
ax.legend()
fig.tight_layout()

import os
plt.savefig(f"python/build/{os.path.basename(__file__).split('.')[0]}.pdf")
plt.show()
