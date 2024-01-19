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

from scipy.optimize import curve_fit

Ts = np.array([0.])
Us = np.array([-2.5])
Vs = np.array([ #"-0.000001", "-0.0000013", "-0.0000015", "-0.0000017", "-0.000002", "-0.0000025", 
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
fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(12.8, 8.2), sharex="col", sharey="row", gridspec_kw=dict(hspace=0, wspace=0))

plotters = np.empty((2,2), dtype=ps.CURVEFAMILY)


for i in range(2):
    counter = 0
    plotters[0][i] = ps.CURVEFAMILY(3, axis=axs[0][i])
    plotters[0][i].set_individual_colors("nice")
    plotters[0][i].set_individual_linestyles(["", "", "-"])
    plotters[0][i].set_individual_markerstyles(["x", "v", ""])
    
    plotters[1][i] = ps.CURVEFAMILY(2, axis=axs[1][i])
    plotters[1][i].set_individual_colors("nice")
    plotters[1][i].set_individual_linestyles(["", "-"])
    plotters[1][i].set_individual_markerstyles(["x", ""])
    
    if i == 1:
        Vs = np.append(Vs, ["-0.3", "-0.34"])
    v_data = np.log(np.array([np.abs(float(v)) for v in Vs]))
    weights = np.zeros(len(Vs))
    peak_positions = np.zeros(len(Vs))
    peak_positions_rel = np.zeros(len(Vs))
    peak_positions_div_delta = np.zeros(len(Vs))
    
    for T, U, V in iterate_containers(Ts, Us, Vs):
        name = f"T={T}/U={U}/V={V}"
        peak_positions[counter], weights[counter] = rp.analyze_peak(f"{folders[i]}{name}", name_suffix)
        peak_positions_div_delta[counter] = peak_positions[counter] / extract_key(f"{folders[i]}{name}/resolvent_{name_suffix}.dat.gz", "Total Gap")
        peak_positions_rel[counter] = 2 * extract_key(f"{folders[i]}{name}/resolvent_{name_suffix}.dat.gz", "Total Gap") - peak_positions[counter]
        counter += 1

    plotters[1][i].plot(v_data, np.exp(weights), label="Data")
    plotters[0][i].plot(v_data, (peak_positions), label=r"$\omega_0$")
    plotters[0][i].plot(v_data, (peak_positions_rel), label=r"$\omega_- - \omega_0$")
    #plotters[0][i].plot(v_data, (peak_positions_div_delta), label=r"$\omega_0 / \Delta_\mathrm{SC}$")
    
    def func(x, a, b, c, d):
        return a * (np.tanh(b * (x + c)) + d)
    
    popt, pcov = curve_fit(func, v_data, np.exp(weights), p0=(4.5, -0.1, 6, 1))
    print(popt)
    lin = np.linspace(min(v_data), max(v_data))
    plotters[1][i].plot(lin, func(lin, *popt), label="Fit")
    
    # V = 0, Delta_CDW = 0
    peak_pos_value, peak_weight = rp.analyze_peak(f"{folders[i][:-1]}_SC/T={Ts[0]}/U={Us[0]}/V=0.0", name_suffix)
    axs[0][i].axhline((peak_pos_value), color="black", linestyle="--")#, label="Value for $V=0$")
    #peak_positions_div_delta[counter] = np.copy(peak_pos_value) / extract_key(f"{folders[i][:-1]}_SC/{name}/resolvent_{name_suffix}.dat.gz", "Total Gap")
    axs[1][i].axhline(np.exp(peak_weight), color="black", linestyle="--")#, label="Value for $V=0$")

axs[0][0].title.set_text("Square")
axs[0][1].title.set_text("Simple cubic")
#axs[0][0].set_ylim(axs[0][0].get_ylim()[0], axs[0][0].get_ylim()[1] + 0.17)
#axs[0][1].set_ylim(axs[0][1].get_ylim()[0], axs[0][1].get_ylim()[1] + 0.17)
#axs[1][0].set_ylim(axs[1][0].get_ylim()[0], axs[1][0].get_ylim()[1] + 0.17)
#axs[1][1].set_ylim(axs[1][1].get_ylim()[0], axs[1][1].get_ylim()[1] + 0.17)

axs[0][0].text(0.88, 0.88, "(a.1)", transform = axs[0][0].transAxes)
axs[0][1].text(0.88, 0.88, "(a.2)", transform = axs[0][1].transAxes)
axs[1][0].text(0.88, 0.88, "(b.1)", transform = axs[1][0].transAxes)
axs[1][1].text(0.88, 0.88, "(b.2)", transform = axs[1][1].transAxes)

axs[1][0].set_xlabel(r"$\ln|V / t|$")
axs[1][1].set_xlabel(r"$\ln|V / t|$")
axs[1][0].set_ylabel(r"$w_0$")
axs[0][0].set_ylabel(r"$\omega_0$")
axs[1][1].legend(loc="lower left")
axs[0][1].legend(loc="upper left")
fig.tight_layout()

plt.savefig(f"plots/{os.path.basename(__file__).split('.')[0]}.pdf")
#plt.show()
