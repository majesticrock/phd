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
                #"-0.00001", "-0.000013", "-0.000015", "-0.000017", "-0.00002", "-0.000025", "-0.00003", "-0.00004", 
                "-0.00005", "-0.00006", "-0.00007", "-0.00008", "-0.00009",  # last index = 25
                "-0.0001", "-0.00013", "-0.00015", "-0.00017", "-0.0002", "-0.00025", "-0.0003", # last index = 32
                "-0.0004", "-0.0005", "-0.0006", "-0.0007", "-0.0008", "-0.0009",  # last index = 37
                "-0.001", "-0.0013", "-0.0015", "-0.0017", "-0.002", "-0.0025", # 43
                "-0.003", "-0.004", "-0.005", "-0.006", "-0.007", "-0.008", "-0.009", # 50
                "-0.01", "-0.013", "-0.015", "-0.017", "-0.02", "-0.025", 
                "-0.03", "-0.04", "-0.05", "-0.06", "-0.07", "-0.08", "-0.09", 
                "-0.1", "-0.13", "-0.15", "-0.2", "-0.25", "-0.28"
                ])

folders = ["../data/modes/square/dos_6000/", "../data/modes/cube/dos_6000/"]
element_names = ["a", "a+b", "a+ib"]

name_suffix = "higgs_CDW"
ncols=2
nrows=2

fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(12.8, 8.2), sharex="col", sharey="row", gridspec_kw=dict(hspace=0, wspace=0))

plotters = np.empty((nrows,ncols), dtype=ps.CURVEFAMILY)
peak_fit_params = {"range": 1e-4, "begin_offset": 1e-10, "imaginary_offset": 1e-7, "peak_position_tol": 1e-14}

cutoffs = [[5, 12], [10, 15]]

for j in range(nrows):
    for i in range(ncols):
        plotters[j][i] = ps.CURVEFAMILY(3, axis=axs[j][i])
        plotters[j][i].set_individual_colors("nice")
        plotters[j][i].set_individual_linestyles(["-", "", ""])
        plotters[j][i].set_individual_markerstyles(["", "^", "v"])

for i in range(2):
    #if i == 1:
    #    Vs = np.append(Vs, ["-0.3", "-0.34"])
    v_data = np.array([float(v) for v in Vs])
    weights = np.zeros(len(Vs))
    peak_positions = np.zeros(len(Vs))
    
    counter = 0
    for T, U, V in iterate_containers(Ts, Us, Vs):
        name = f"T={T}/U={U}/V={V}"
        peak_positions[counter], weights[counter] = rp.find_weight(f"{folders[i]}{name}", name_suffix, **peak_fit_params)
        counter += 1
    
    v_data = np.log(np.abs(v_data))
    peak_positions = np.log(peak_positions)
    # Plot and fit omega_0
    cut = len(v_data) - cutoffs[0][i]
    popt, pcov = ez_linear_fit(v_data[:cut], peak_positions[:cut], plotters[0][i], ez_lin_space(v_data), label="Fit")
    plotters[0][i].plot(v_data[:cut], peak_positions[:cut], label="Fitted data")
    plotters[0][i].plot(v_data[cut:], peak_positions[cut:], label="Omitted data", markerfacecolor="None")
    axs[0][i].text(0.15, 0.88, f"$a={popt[0]:.4f}$", transform = axs[0][i].transAxes)
    axs[0][i].text(0.15, 0.77, f"$b={popt[1]:.4f}$", transform = axs[0][i].transAxes)
    
    # Plot and fit W_0
    cut = len(v_data) - cutoffs[1][i]
    popt, pcov = ez_linear_fit(v_data[:cut], weights[:cut], plotters[1][i], ez_lin_space(v_data), label="Fit")
    plotters[1][i].plot(v_data[:cut], weights[:cut], label="Data")
    plotters[1][i].plot(v_data[cut:], weights[cut:], label="Omitted data", markerfacecolor="None")
    axs[1][i].text(0.15, 0.88, f"$a={popt[0]:.4f}$", transform = axs[1][i].transAxes)
    axs[1][i].text(0.15, 0.77, f"$b={popt[1]:.4f}$", transform = axs[1][i].transAxes)

axs[0][0].title.set_text("Square")
axs[0][1].title.set_text("Simple cubic")

for i in range(ncols):
    axs[0][i].text(0.88, 0.75, f"(a.{i+1})", transform = axs[0][i].transAxes)
    axs[1][i].text(0.88, 0.75, f"(b.{i+1})", transform = axs[1][i].transAxes)

axs[1][0].set_ylim(axs[1][0].get_ylim()[0], axs[1][0].get_ylim()[1] + 1.1)

axs[1][0].set_xlabel(r"$\ln|V / t|$")
axs[1][1].set_xlabel(r"$\ln|V / t|$")

axs[0][0].set_ylabel(r"$\ln (\omega_0 / t) $")
axs[1][0].set_ylabel(r"$\ln(W_0)$")

#axs[0][1].legend(loc="lower right")
axs[1][1].legend(loc="lower left")
fig.tight_layout()

plt.savefig(f"plots/{os.path.basename(__file__).split('.')[0]}.pdf")
#plt.show()
