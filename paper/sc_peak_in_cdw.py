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
Vs = np.array([#"0.00001", "0.000013", "0.000015", "0.000017", "0.00002", "0.000025", "0.00003", "0.00004", 
                "0.00005", "0.00006", "0.00007", "0.00008", "0.00009", 
                
                "0.0001", "0.00013", "0.00015", "0.00017", "0.0002", "0.00025", "0.0003", 
                "0.0004", "0.0005", "0.0006", "0.0007", "0.0008", "0.0009", 
                
                "0.001", "0.0013", "0.0015", "0.0017", "0.002", "0.0025", 
                "0.003", "0.004", "0.005", "0.006", "0.007", "0.008", "0.009", 
                
                "0.01", "0.013", "0.015", "0.017", "0.02", "0.025", 
                "0.03", "0.04", "0.05", "0.06", "0.07", "0.08", "0.09", 
                
                "0.1", "0.13", "0.15", "0.2", "0.25", "0.3", "0.35", "0.4", 
                "0.45", "0.5", "0.6", "0.7", "0.8", "0.9", 
                "1.0", "1.1", "1.2", "1.3", "1.4", "1.5", "1.6", "1.7", "1.8", "1.9", "2.0", 
                "2.5", "3.0", "3.5", "4.0", "6.0", "8.0", "10.0", "15.0", "25.0", "50.0"])
v_data = np.log(np.array([float(v) for v in Vs]))

folders = ["../data/modes/square/dos_6000/", "../data/modes/cube/dos_6000/"]
element_names = ["a", "a+b", "a+ib"]

name_suffix = "phase_SC"

weights = np.zeros(len(Vs))
peak_positions = np.zeros(len(Vs))
peak_positions_div_delta = np.zeros(len(Vs))
fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(12.8, 8.2), sharex=True, sharey="row", gridspec_kw=dict(hspace=0, wspace=0))

plotters = np.empty((2,2), dtype=ps.CURVEFAMILY)


for i in range(2):
    counter = 0
    plotters[0][i] = ps.CURVEFAMILY(6, axis=axs[0][i])
    plotters[0][i].set_individual_colors("nice")
    plotters[0][i].set_individual_linestyles(["-", "-", "", "", "", ""])
    plotters[0][i].set_individual_markerstyles(["", "", "^", "^", "v", "1"])
    
    plotters[1][i] = ps.CURVEFAMILY(3, axis=axs[1][i])
    plotters[1][i].set_individual_colors("nice")
    plotters[1][i].set_individual_linestyles(["-", "", ""])
    plotters[1][i].set_individual_markerstyles(["", "^", "v"])
    
    for T, U, V in iterate_containers(Ts, Us, Vs):
        name = f"T={T}/U={U}/V={V}"
        cont_edges = cf.continuum_edges(f"{folders[i]}{name}", name_suffix, True)
        lower = 0 if float(V) < 1 else float(V)
        upper = 13 * float(V) + 2 if 13 * float(V) + 2 < cont_edges[0] else cont_edges[0]
        
        peak_positions[counter], weights[counter] = rp.analyze_peak(f"{folders[i]}{name}", name_suffix, (lower, upper), begin_offset=1e-7)
        peak_positions_div_delta[counter] = peak_positions[counter] / extract_key(f"{folders[i]}{name}/resolvent_{name_suffix}.dat.gz", "Total Gap")
        counter += 1

    cut = 38
    popt, pcov = curve_fit(rp.linear_function, v_data[:cut], weights[:cut])
    v_lin = np.linspace(v_data.min(), v_data.max(), 500)
    
    plotters[1][i].plot(v_lin, rp.linear_function(v_lin, *popt), label="Fit")
    plotters[1][i].plot(v_data[:cut], weights[:cut], label="Fitted data")
    plotters[1][i].plot(v_data[cut:], weights[cut:], label="Omitted data", markerfacecolor="None")
    
    axs[1][i].text(0.05, 0.4, f"$c={popt[0]:.4f}$", transform = axs[1][i].transAxes)
    axs[1][i].text(0.05, 0.3, f"$d={popt[1]:.4f}$", transform = axs[1][i].transAxes)
    
    peak_positions = np.log(peak_positions)
    peak_positions_div_delta = np.log(peak_positions_div_delta)
    popt, pcov = curve_fit(rp.linear_function, v_data[:cut], peak_positions[:cut])
    x_lin = np.linspace(np.min(v_data), np.max(v_data), 2000)
    plotters[0][i].plot(x_lin, rp.linear_function(x_lin, *popt), label="Fit 1")
    axs[0][i].text(0.05, 0.9, f"$c_1={popt[0]:.4f}$", transform = axs[0][i].transAxes)
    axs[0][i].text(0.05, 0.8, f"$d_1={popt[1]:.4f}$", transform = axs[0][i].transAxes)

    cut2 = len(Vs) - 27
    popt, pcov = curve_fit(rp.linear_function, v_data[cut2:], peak_positions[cut2:])
    axs[0][i].text(0.05, 0.7, f"$c_2={popt[0]:.4f}$", transform = axs[0][i].transAxes)
    axs[0][i].text(0.05, 0.6, f"$d_2={popt[1]:.4f}$", transform = axs[0][i].transAxes)
    plotters[0][i].plot(x_lin,  rp.linear_function(x_lin, *popt), label="Fit 2")
    plotters[0][i].plot(v_data[:cut], peak_positions[:cut], label="Data Fit 1")
    plotters[0][i].plot(v_data[cut2:], peak_positions[cut2:], label="Data Fit 2")
    plotters[0][i].plot(v_data[cut:cut2], peak_positions[cut:cut2], label="Omitted data", markerfacecolor="None")
    plotters[0][i].plot(v_data, peak_positions_div_delta, label=r"$\omega_0 / \Delta_\mathrm{CDW}$")

axs[0][0].title.set_text("Square")
axs[0][1].title.set_text("Simple cubic")

axs[0][0].text(0.8, 0.92, "(a.1)", transform = axs[0][0].transAxes)
axs[0][1].text(0.8, 0.92, "(a.2)", transform = axs[0][1].transAxes)
axs[1][0].text(0.8, 0.92, "(b.1)", transform = axs[1][0].transAxes)
axs[1][1].text(0.8, 0.92, "(b.2)", transform = axs[1][1].transAxes)

axs[1][0].set_xlabel(r"$\ln(V / t)$")
axs[1][1].set_xlabel(r"$\ln(V / t)$")
axs[1][0].set_ylabel(r"$b = \ln(W_0)$")
axs[0][0].set_ylabel(r"$\ln(\omega_0 / t)$")
axs[1][1].legend(loc="lower right")
axs[0][1].legend(loc="lower right", ncol=2)
fig.tight_layout()

plt.savefig(f"plots/{os.path.basename(__file__).split('.')[0]}.pdf")
#plt.show()
