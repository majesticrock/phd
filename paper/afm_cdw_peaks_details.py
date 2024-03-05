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
Us_square = np.array([  0.0001, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 
                        0.8, 0.85, 0.9, 0.95, 1, 1.05, 1.1, 1.15, 1.2, 1.24, 1.25, 1.27, 1.3, 1.32, 1.35, 1.37, 1.4, 1.42])
#new: 1.422, 1.424, 1.426, 1.428, 1.43, 1.432, 1.434, 1.436, 1.438, 1.44
Us_cube = np.array([0.0001, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 
                    0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.295, 0.299, 0.3, 
                    0.302, 0.305, 0.307, 0.31, 0.312, 0.315, 0.317, 0.32])#~0.3201
#new 0.321

Us = np.array([Us_square, Us_cube], dtype=object)
Vs = np.array([0.])


folders = ["../data/modes/square/dos_6000/", "../data/modes/cube/dos_6000/"]
element_names = ["a", "a+b", "a+ib"]

name_suffix = "AFM"
nrows=2
fig, axs = plt.subplots(nrows=nrows, ncols=2, figsize=(12.8, 5.4), sharex="col", gridspec_kw=dict(hspace=0))
plotters = np.empty((nrows,2), dtype=ps.CURVEFAMILY)

peak_fit_params = [ {"range": 1e-6, "begin_offset": 5e-10, "imaginary_offset": 1e-6},
                    {"range": 5e-8, "begin_offset": 1e-12, "imaginary_offset": 1e-8}]

for i in range(2):
    for j in range(nrows):
        plotters[j][i] = ps.CURVEFAMILY(3, axis=axs[j][i])
        plotters[j][i].set_individual_colors("nice")
        plotters[j][i].set_individual_linestyles(["-", "", ""])
        plotters[j][i].set_individual_markerstyles(["", "^", "v"])
    u_data = 4.8 - Us[i]
    
    weights = np.zeros(len(u_data))
    peak_positions = np.zeros(len(u_data))
    peak_positions_to_cont = np.zeros(len(u_data))
    counter = 0
    for T, U, V in iterate_containers(Ts, u_data, Vs):
        V = 1.2 if i == 0 else 0.8
        name = f"T={T}/U={round(U, 4)}/V={V}"
        cont_edges = cf.continuum_edges(f"{folders[i]}{name}", f"higgs_{name_suffix}", True)
        lower = 0 if float(V) < 1 else float(V)
        upper = 13 * float(V) + 2 if 13 * float(V) + 2 < cont_edges[0] else cont_edges[0]

        peak_positions[counter], weights[counter] = rp.find_weight(f"{folders[i]}{name}", f"higgs_{name_suffix}", (lower, upper), peak_position_tol=1e-14, reversed=True, **(peak_fit_params[0 if counter < len(u_data) - 5 else 1]))
        peak_positions_to_cont[counter] = cont_edges[0] - peak_positions[counter]
        counter += 1
    
    u_log = np.log(u_data - min(u_data) + 0.008)
    peak_positions_to_cont = np.log(peak_positions_to_cont)
    
    cut = 15
    popt1, pcov1 = ez_linear_fit(u_log[cut:], peak_positions_to_cont[cut:], plotters[0][i], ez_lin_space(u_log), label="Fit")
    popt2, pcov2 = ez_linear_fit(u_log[cut:], weights[cut:], plotters[1][i], ez_lin_space(u_log), label="Fit")  
    plotters[0][i].plot(u_log[cut:], peak_positions_to_cont[cut:], label=f"Fitted data")
    plotters[1][i].plot(u_log[cut:], weights[cut:])
    
    plotters[0][i].plot(u_log[:cut], peak_positions_to_cont[:cut], label=f"Omitted data", markerfacecolor="None")
    plotters[1][i].plot(u_log[:cut], weights[:cut], markerfacecolor="None")
    
    axs[0][i].text(0.05, 0.8, f"a = {popt1[0]:.4f}", transform = axs[0][i].transAxes)
    axs[1][i].text(0.05, 0.8, f"a = {popt2[0]:.4f}", transform = axs[1][i].transAxes)
    axs[0][i].text(0.05, 0.7, f"b = {popt1[1]:.4f}", transform = axs[0][i].transAxes)
    axs[1][i].text(0.05, 0.7, f"b = {popt2[1]:.4f}", transform = axs[1][i].transAxes)

axs[0][0].title.set_text("Square lattice")
axs[0][1].title.set_text("Simple cubic lattice")

axs[0][0].text(0.89, 0.7, "(a.1)", transform = axs[0][0].transAxes)
axs[0][1].text(0.89, 0.7, "(a.2)", transform = axs[0][1].transAxes)
axs[1][0].text(0.89, 0.7, "(b.1)", transform = axs[1][0].transAxes)
axs[1][1].text(0.89, 0.7, "(b.2)", transform = axs[1][1].transAxes)

axs[nrows-1][0].set_xlabel(r"$\ln((U - U_0) / t)$")
axs[nrows-1][1].set_xlabel(r"$\ln((U - U_0) / t)$")
axs[0][0].set_ylabel(r"$\ln ((\omega_- - \omega_0) / t)$")
axs[1][0].set_ylabel(r"$\ln(W_0)$")
legend = axs[0][1].legend(loc='lower right')
fig.tight_layout()

plt.savefig(f"plots/{os.path.basename(__file__).split('.')[0]}.pdf")
#plt.show()
