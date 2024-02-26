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


Ts = np.array([0.])
Us = np.array([-2.5])
Vs = np.array([ #"-0.000001", "-0.0000013", "-0.0000015", "-0.0000017", "-0.000002", "-0.0000025", 
                #"-0.000003", "-0.000004", "-0.000005", "-0.000006", "-0.000007", "-0.000008", "-0.000009", # last index = 12
                "-0.00001", "-0.000013", "-0.000015", "-0.000017", "-0.00002", "-0.000025", 
                "-0.00003", "-0.00004", "-0.00005", "-0.00006", "-0.00007", "-0.00008", "-0.00009",  # last index = 25
                "-0.0001", "-0.00013", "-0.00015", "-0.00017", "-0.0002", "-0.00025", "-0.0003", # last index = 32
                "-0.0004", "-0.0005", "-0.0006", "-0.0007", "-0.0008", "-0.0009",  # last index = 37
                "-0.001", "-0.0013", "-0.0015", "-0.0017", "-0.002", "-0.0025", # 43
                "-0.003", "-0.004", "-0.005", "-0.006", "-0.007", "-0.008", "-0.009", # 50
                "-0.01", "-0.013", "-0.015", "-0.017", "-0.02", "-0.025", 
                "-0.03", "-0.04", "-0.05", "-0.06", "-0.07", "-0.08", "-0.09", 
                "-0.1", "-0.13", "-0.15", "-0.2", "-0.25", "-0.28"
                ])

# Data points not to plot for the non-log plot
omit_non_log = [1, 
                #2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 
                #13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
                #27, 28, 29, 30, 31, 32, 33, 
                #35, 36, 37, 
                #39, 40, 41, 42, 43,
                #45, 47, 48, 50, 52, 54
                ]

folders = ["../data/modes/square/dos_6k/", "../data/modes/cube/dos_3k/"]
element_names = ["a", "a+b", "a+ib"]

name_suffix = "higgs_CDW"
ncols=2
nrows=4

fig = plt.figure(figsize=(12.8, 10))
gs_outer = fig.add_gridspec(2, 1)
gs_top = gs_outer[0, 0].subgridspec(2, 2, wspace=0, hspace=0)
gs_bot = gs_outer[1, 0].subgridspec(2, 2, wspace=0, hspace=0)

axs_top = gs_top.subplots(sharey="row", sharex="col")
axs_bot = gs_bot.subplots(sharey="row", sharex="col")
axs = np.vstack((axs_top, axs_bot))

for j in range(ncols):
    axs[0][j].sharex(axs[1][j])
    axs[2][j].sharex(axs[3][j])

plotters = np.empty((nrows,ncols), dtype=ps.CURVEFAMILY)
peak_fit_params = {"range": 1e-4, "begin_offset": 1e-10, "imaginary_offset": 1e-7, "peak_position_tol": 1e-14}

for j in range(nrows):
    for i in range(ncols):
        plotters[j][i] = ps.CURVEFAMILY(2, axis=axs[j][i])
        plotters[j][i].set_individual_colors("nice")
        plotters[j][i].set_individual_linestyles(["-", ""])
        plotters[j][i].set_individual_markerstyles(["", "^"])

for i in range(2):
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
    
    ## V = 0, Delta_CDW = 0
    #peak_pos_value, peak_weight = rp.find_weight(f"{folders[i][:-1]}_SC/T={Ts[0]}/U={Us[0]}/V=0.0", name_suffix, **peak_fit_params)
    original_data = np.array([peak_positions, np.exp(weights)])
    original_v = np.copy(v_data)
    #
    #peak_positions = (peak_positions - peak_pos_value)
    #weights = (np.exp(peak_weight) - np.exp(weights))
    
    def test_func(x, a, b, c):
        return a * np.tanh(b * x - c) + a
    v_data = np.log(np.abs(v_data))
    
    # Plot and fit omega_0
    #popt, pcov = ez_general_fit(v_data, peak_positions, test_func, plotters[2][i], ez_lin_space(v_data), label="Fit")
    plotters[2][i].plot(v_data, np.log(peak_positions), label="Data")
    #axs[2][i].text(0.05, 0.88, f"$a={popt[0]:.4f}$", transform = axs[2][i].transAxes)
    #axs[2][i].text(0.05, 0.77, f"$b={popt[1]:.4f}$", transform = axs[2][i].transAxes)
    #axs[2][i].text(0.05, 0.66, f"$c={popt[2]:.4f}$", transform = axs[2][i].transAxes)
    
    v_lin = ez_lin_space(original_v)
    v_log = np.log(np.abs(v_lin))
    original_v = np.delete(original_v, omit_non_log)
    
    #plotters[0][i].plot(v_lin, test_func(v_log, *popt) + peak_pos_value, label="Fit")
    plotters[0][i].plot(original_v, np.delete(original_data[0], omit_non_log), label="Data")
    #axs[0][i].axhline(peak_pos_value, color="k", ls="--")
    
    # Plot and fit W_0
    #popt, pcov = ez_general_fit(v_data, weights, test_func, plotters[3][i], ez_lin_space(v_data), label="Fit")
    plotters[3][i].plot(v_data, weights, label="Data")
    #axs[3][i].text(0.05, 0.88, f"$a={popt[0]:.4f}$", transform = axs[3][i].transAxes)
    #axs[3][i].text(0.05, 0.77, f"$b={popt[1]:.4f}$", transform = axs[3][i].transAxes)
    #axs[3][i].text(0.05, 0.66, f"$c={popt[2]:.4f}$", transform = axs[3][i].transAxes)
    
    #plotters[1][i].plot(v_lin, -test_func(v_log, *popt) + np.exp(peak_weight), label="Fit")
    plotters[1][i].plot(original_v, np.delete(original_data[1], omit_non_log), label="Data")
    #axs[1][i].axhline(np.exp(peak_weight), color="k", ls="--")

axs[0][0].title.set_text("Square")
axs[0][1].title.set_text("Simple cubic")

for i in range(ncols):
    axs[0][i].text(0.88, 0.85 - 0.09*i, f"(a.{i+1})", transform = axs[0][i].transAxes)
    axs[1][i].text(0.88, 0.85 - 0.09*i, f"(b.{i+1})", transform = axs[1][i].transAxes)
    axs[2][i].text(0.88, 0.85 - 0.06*i, f"(c.{i+1})", transform = axs[2][i].transAxes)
    axs[3][i].text(0.88, 0.85 - 0.06*i, f"(d.{i+1})", transform = axs[3][i].transAxes)

axs[1][0].set_ylim(axs[1][0].get_ylim()[0], axs[1][0].get_ylim()[1] + 1.1)

axs[1][0].set_xlabel(r"$V [t]$")
axs[1][1].set_xlabel(r"$V [t]$")
axs[3][0].set_xlabel(r"$\ln|V / t|$")
axs[3][1].set_xlabel(r"$\ln|V / t|$")

axs[0][0].set_ylabel(r"$\omega_0 [t]$")
axs[1][0].set_ylabel(r"$W_0$")
axs[2][0].set_ylabel(r"$(\omega_0 - \omega_0^{V = 0}) [t]$")
axs[3][0].set_ylabel(r"$W_0^{V = 0} - W_0$")

#axs[0][1].legend(loc="lower right")
axs[1][1].legend(loc="upper left")
fig.tight_layout()

plt.savefig(f"plots/{os.path.basename(__file__).split('.')[0]}.pdf")
#plt.show()
