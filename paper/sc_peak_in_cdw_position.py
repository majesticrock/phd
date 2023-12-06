import numpy as np
import matplotlib.pyplot as plt

import os, sys
if os.name == "nt":
    sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + r"\python")
else:
    sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + "/python")

from lib.iterate_containers import *
from lib.extract_key import *
import lib.resolvent_peak as rp
import lib.plot_settings as ps

from scipy.optimize import curve_fit

Ts = np.array([0.])
Us = np.array([-2.5])
Vs = np.array(["50.0", "25.0", "15.0", "10.0", "8.0", "6.0", "4.0", "3.5", "3.0", "2.5", "2.0", "1.5", "1.4", "1.3",
               "1.2", "1.1", "1.0", "0.9", "0.8", "0.7", "0.6", "0.5", 
               "0.45", "0.4", "0.35", "0.3", "0.25", "0.2", "0.15",
               "0.1", "0.07", "0.05", "0.04", "0.03", "0.02", "0.01", 
               "0.007", "0.005", "0.004", "0.003", "0.002", "0.0015", "0.001", 
               "0.0007", "0.0005", "0.0004", "0.0003", "0.0002", "0.00015", "0.0001", 
               "0.00007", "0.00005", "0.00003"])
v_data = np.log(np.array([float(v) for v in Vs]))

folders = ["../data/modes/square/dos_3k/", "../data/modes/cube/dos_3k/"]
element_names = ["a", "a+b", "a+ib"]

name_suffix = "phase_SC"

weights = np.zeros(len(Vs))
peak_positions = np.zeros(len(Vs))
fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(12.8, 8.2), sharex=True, sharey="row", gridspec_kw=dict(hspace=0, wspace=0))

plotters = np.empty((2,2), dtype=ps.CURVEFAMILY)

for i in range(2):
    counter = 0
    plotters[0][i] = ps.CURVEFAMILY(5, axis=axs[0][i])
    plotters[0][i].set_individual_colors("nice")
    plotters[0][i].set_individual_linestyles(["-", "-", "", "", ""])
    plotters[0][i].set_individual_markerstyles(["", "", "X", "X", "o"])
    
    plotters[1][i] = ps.CURVEFAMILY(3, axis=axs[1][i])
    plotters[1][i].set_individual_colors("nice")
    plotters[1][i].set_individual_linestyles(["-", "", ""])
    plotters[1][i].set_individual_markerstyles(["", "X", "o"])
    
    for T, U, V in iterate_containers(Ts, Us, Vs):
        name = f"T={T}/U={U}/V={V}"
        lower = float(V)
        upper = 8 * float(V) + 2

        peak = rp.Peak(f"{folders[i]}{name}", name_suffix, (lower, upper))
        peak_pos_value = np.copy(peak.peak_position)
        scipy_result = peak.improved_peak_position(1e-2, 1e-12)
        # only an issue if the difference is too large;
        if scipy_result[2]["warnflag"] != 0 and np.abs((scipy_result[0][0] - peak_pos_value) / peak_pos_value) > 1e-3:
            print(f"We might not have found the peak for V={V}!\nWe found ", peak_pos_value, " and\n", scipy_result)

        peak_pos_value = scipy_result[0][0]
        peak_positions[counter] = np.copy(scipy_result[0][0])

        popt, pcov, w_log, y_data = peak.fit_real_part(0.01, 1e-8)
        weights[counter] = popt[1]
        counter += 1

    cut = -18
    popt, pcov = curve_fit(rp.linear_function, v_data[cut:], weights[cut:])
    v_lin = np.linspace(v_data.min(), v_data.max(), 500)
    
    plotters[1][i].plot(v_lin, rp.linear_function(v_lin, *popt), label="Fit")
    plotters[1][i].plot(v_data[cut:], weights[cut:], "X", label="Fitted data")
    plotters[1][i].plot(v_data[:cut], weights[:cut], "o", label="Omitted data")
    
    axs[1][i].text(0.05, 0.4, f"$c={popt[0]:.4f}\pm{np.sqrt(pcov[0][0]):.4f}$", transform = axs[1][i].transAxes)
    axs[1][i].text(0.05, 0.3, f"$d={popt[1]:.4f}\pm{np.sqrt(pcov[1][1]):.4f}$", transform = axs[1][i].transAxes)
    
    peak_positions = np.log(peak_positions)
    popt, pcov = curve_fit(rp.linear_function, v_data[cut:], peak_positions[cut:])
    x_lin = np.linspace(np.min(v_data), np.max(v_data), 2000)
    plotters[0][i].plot(x_lin, rp.linear_function(x_lin, *popt), label="Fit 1")
    axs[0][i].text(0.6, 0.4, f"$c_1={popt[0]:.4f}\pm{np.sqrt(pcov[0][0]):.4f}$", transform = axs[0][i].transAxes)
    axs[0][i].text(0.6, 0.3, f"$d_1={popt[1]:.4f}\pm{np.sqrt(pcov[1][1]):.4f}$", transform = axs[0][i].transAxes)

    cut2 = 20
    popt, pcov = curve_fit(rp.linear_function, v_data[:cut2], peak_positions[:cut2])
    axs[0][i].text(0.6, 0.2, f"$c_2={popt[0]:.4f}\pm{np.sqrt(pcov[0][0]):.4f}$", transform = axs[0][i].transAxes)
    axs[0][i].text(0.6, 0.1, f"$d_2={popt[1]:.4f}\pm{np.sqrt(pcov[1][1]):.4f}$", transform = axs[0][i].transAxes)
    plotters[0][i].plot(x_lin,  rp.linear_function(x_lin, *popt), label="Fit 2")
    plotters[0][i].plot(v_data[cut:], peak_positions[cut:], "X", label="Data Fit 1")
    plotters[0][i].plot(v_data[:cut2], peak_positions[:cut2], "X", label="Data Fit 2")
    plotters[0][i].plot(v_data[cut2:cut], peak_positions[cut2:cut], "o", label="Omitted data")

axs[0][0].title.set_text("Square")
axs[0][1].title.set_text("Simple cubic")

axs[1][0].set_xlabel(r"$\ln(V / t)$")
axs[1][1].set_xlabel(r"$\ln(V / t)$")
axs[1][0].set_ylabel(r"$b = \ln(w_0 \cdot t)$")
axs[0][0].set_ylabel(r"$\ln(\omega_0 / t)$")
axs[1][0].legend(loc="lower right")
axs[0][0].legend(loc="upper left")
fig.tight_layout()

#axs[1][0].grid()
#axs[1][1].grid()
#axs[0][0].grid()
#axs[0][1].grid()

import os
plt.savefig(f"plots/{os.path.basename(__file__).split('.')[0]}.pdf")
#plt.show()
