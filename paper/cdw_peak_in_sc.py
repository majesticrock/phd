import numpy as np
import matplotlib.pyplot as plt

import os, sys
if os.name == "nt":
    sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + r"\python")
else:
    sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + "/python")

import lib.continued_fraction as cf
from lib.iterate_containers import *
from lib.extract_key import *
import lib.resolvent_peak as rp
import lib.plot_settings as ps
from lib.ez_fit import *

from scipy.optimize import curve_fit

Ts = np.array([0.])
Us = np.array([-2.5])
Vs = np.array([ "-0.000001", "-0.0000013", "-0.0000015", "-0.0000017", "-0.000002", "-0.0000025", 
                "-0.000003", "-0.000004", "-0.000005", "-0.000006", "-0.000007", "-0.000008", "-0.000009",
                "-0.00001", "-0.000013", "-0.000015", "-0.000017", "-0.00002", "-0.000025", 
                "-0.00003", "-0.00004", "-0.00005", "-0.00006", "-0.00007", "-0.00008", "-0.00009", 
                "-0.0001", "-0.00013", "-0.00015", "-0.00017", "-0.0002", "-0.00025", "-0.0003", 
                "-0.0004", "-0.0005", "-0.0006", "-0.0007", "-0.0008", "-0.0009", 
                "-0.001", "-0.0013", "-0.0015", "-0.0017", "-0.002", "-0.0025", 
                "-0.003", "-0.004", "-0.005", "-0.006", "-0.007", "-0.008", "-0.009", 
                "-0.01", "-0.013", "-0.015", "-0.017", "-0.02", "-0.025", 
                "-0.03", "-0.04", "-0.05", "-0.06", "-0.07", "-0.08", "-0.09", 
                "-0.1", "-0.13", "-0.15", "-0.2", "-0.25", "-0.28"
                ])

folders = ["../data/modes/square/dos_3k/", "../data/modes/cube/dos_3k/"]
element_names = ["a", "a+b", "a+ib"]

name_suffix = "higgs_CDW"
fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(12.8, 6.4), sharex="col", sharey="row", gridspec_kw=dict(hspace=0, wspace=0))

plotters = np.empty((2,2), dtype=ps.CURVEFAMILY)
peak_fit_params = {"range": 1e-4, "begin_offset": 1e-10, "imaginary_offset": 1e-7, "peak_position_tol": 1e-14}

for i in range(2):
    for j in range(2):
        plotters[j][i] = ps.CURVEFAMILY(2, axis=axs[j][i])
        plotters[j][i].set_individual_colors("nice")
        plotters[j][i].set_individual_linestyles(["-", ""])
        plotters[j][i].set_individual_markerstyles(["", "x"])

    if i == 1:
        Vs = np.append(Vs, ["-0.3", "-0.34"])
    v_data = np.array([float(v) for v in Vs])
    weights = np.zeros(len(Vs))
    peak_positions = np.zeros(len(Vs))
    
    counter = 0
    for T, U, V in iterate_containers(Ts, Us, Vs):
        name = f"T={T}/U={U}/V={V}"
        peak_positions[counter], weights[counter] = rp.find_weight(f"{folders[i]}{name}", name_suffix, **peak_fit_params)
        counter += 1
    
    # V = 0, Delta_CDW = 0
    peak_pos_value, peak_weight = rp.find_weight(f"{folders[i][:-1]}_SC/T={Ts[0]}/U={Us[0]}/V=0.0", name_suffix, **peak_fit_params)
    peak_positions = (peak_positions - peak_pos_value)
    weights = (np.exp(peak_weight) - np.exp(weights))
    
    def test_func(x, a, b, c):
        return a * np.tanh(b * x - c) + a
    v_data = np.log(np.abs(v_data))
    
    # Plot and fit omega_0
    popt, pcov = ez_general_fit(v_data, peak_positions, test_func, plotters[0][i], ez_lin_space(v_data), label="Fit")
    plotters[0][i].plot(v_data, peak_positions, label="Fitted data")
    axs[0][i].text(0.05, 0.9, f"$a={popt[0]:.4f}$", transform = axs[0][i].transAxes)
    axs[0][i].text(0.05, 0.8, f"$b={popt[1]:.4f}$", transform = axs[0][i].transAxes)
    axs[0][i].text(0.05, 0.7, f"$c={popt[2]:.4f}$", transform = axs[0][i].transAxes)
    
    # Plot and fit W_0
    popt, pcov = ez_general_fit(v_data, weights, test_func, plotters[1][i], ez_lin_space(v_data), label="Fit")
    plotters[1][i].plot(v_data, weights, label="Fitted data")
    axs[1][i].text(0.05, 0.9, f"$a={popt[0]:.4f}$", transform = axs[1][i].transAxes)
    axs[1][i].text(0.05, 0.8, f"$b={popt[1]:.4f}$", transform = axs[1][i].transAxes)
    axs[1][i].text(0.05, 0.7, f"$c={popt[2]:.4f}$", transform = axs[1][i].transAxes)

axs[0][0].title.set_text("Square")
axs[0][1].title.set_text("Simple cubic")

axs[0][0].text(0.88, 0.88, "(a.1)", transform = axs[0][0].transAxes)
axs[0][1].text(0.88, 0.88, "(a.2)", transform = axs[0][1].transAxes)
axs[1][0].text(0.88, 0.88, "(b.1)", transform = axs[1][0].transAxes)
axs[1][1].text(0.88, 0.88, "(b.2)", transform = axs[1][1].transAxes)

axs[1][0].set_xlabel(r"$\ln|V / t|$")
axs[1][1].set_xlabel(r"$\ln|V / t|$")
axs[0][0].set_ylabel(r"$\ln(\omega_0 - \omega_0(V = 0)) / t$")
axs[1][0].set_ylabel(r"$\ln(W_0(V=0) - W_0)$")
axs[0][1].legend(loc="lower right")
axs[1][1].legend(loc="lower right")
fig.tight_layout()

plt.savefig(f"plots/{os.path.basename(__file__).split('.')[0]}.pdf")
#plt.show()
