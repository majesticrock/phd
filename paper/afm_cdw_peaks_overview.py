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

Us_square = np.array([  0.0001, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 
                        0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.428])#, 1.43, 1.432, 1.434, 1.436, 1.438, 1.44, 1.45, 1.47, 1.5

Us_cube = np.array([0.0001, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 
                    0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31, 0.321])#, 0.322, 0.325

#for p in [0.321]:
#    print(f"0 {round(4.8 + p, 4)} 1.2")

Ts = np.array([0.])
Us = np.array([np.concatenate((-Us_square[::-1], Us_square)), np.concatenate((-Us_cube[::-1], Us_cube))], dtype=object)
Vs = np.array([1.0])

folders = ["../data/modes/square/dos_6000/", "../data/modes/cube/dos_6000/"]
element_names = ["a", "a+b", "a+ib"]

name_suffices = ["CDW", "AFM"]
nrows=3
fig, axs = plt.subplots(nrows=nrows, ncols=2, figsize=(12.8, 6.4), sharex="col", gridspec_kw=dict(hspace=0))
plotters = np.empty((nrows,2), dtype=ps.CURVEFAMILY)

peak_fit_params = [ {"range": 1e-7, "begin_offset": 1e-10, "imaginary_offset": 5e-7},
                    {"range": 1e-8, "begin_offset": 1e-12, "imaginary_offset": 1e-8}]

for i in range(2):
    for j in range(nrows):
        plotters[j][i] = ps.CURVEFAMILY(6, axis=axs[j][i])
        plotters[j][i].set_individual_colors(["C1", "C2"])
        plotters[j][i].set_individual_linestyles(["-", "-"])
        plotters[j][i].set_individual_markerstyles(["v", "^"])

    u_data = 4.8 + Us[i]
    weights = np.zeros(len(u_data))
    peak_positions = np.zeros(len(u_data))
    peak_positions_to_cont = np.zeros(len(u_data))
    
    for name_suffix in name_suffices:
        counter = 0
        for T, U, V_a in iterate_containers(Ts, u_data, Vs):
            V = 1.2 if i == 0 else 0.8
            name = f"T={T}/U={round(U, 4)}/V={V}"
            cont_edges = cf.continuum_edges(f"{folders[i]}{name}", f"higgs_{name_suffix}", True)
            lower = 0 if float(V) < 1 else float(V)
            upper = 13 * float(V) + 2 if 13 * float(V) + 2 < cont_edges[0] else cont_edges[0]

            peak_positions[counter], weights[counter] = rp.analyze_peak(f"{folders[i]}{name}", f"higgs_{name_suffix}", (lower, upper), peak_position_tol=1e-14, 
                                                                        reversed=True, **(peak_fit_params[1 if counter > len(u_data) - 6 or counter < 5 else 0]))
            peak_positions_to_cont[counter] = cont_edges[0] - peak_positions[counter]
            counter += 1
        
        label_subscript = name_suffix if name_suffix != "AFM" else "l.AFM"
        
        #us_same = np.array([np.concatenate((Us_square[::-1], Us_square)), np.concatenate((Us_cube[::-1], Us_cube))], dtype=object)
        plotters[0][i].plot(u_data, peak_positions, label=f"$\\mathcal{{A}}_\\mathrm{{{label_subscript}}} (\\omega)$")
        plotters[1][i].plot(u_data, peak_positions_to_cont)
        plotters[2][i].plot(u_data, np.exp(weights))
    
    for j in range(nrows):
        axs[j][i].axvspan(min(u_data), (min(u_data) + max(u_data)) / 2, alpha=0.3, color="C1")
        axs[j][i].axvspan((min(u_data) + max(u_data)) / 2, max(u_data), alpha=0.3, color="C2")

axs[0][0].set_title("Square lattice")
axs[0][1].set_title("Simple cubic lattice")

for i in range(2):
    axs[0][i].text(0.87, 0.5, f"(a.{i+1})", transform = axs[0][i].transAxes)
    axs[1][i].text(0.87, 0.5, f"(b.{i+1})", transform = axs[1][i].transAxes)
    axs[2][i].text(0.87, 0.5, f"(c.{i+1})", transform = axs[2][i].transAxes)

axs[nrows-1][0].set_xlabel(r"$U [t]$")
axs[nrows-1][1].set_xlabel(r"$U [t]$")
axs[0][0].set_ylabel(r"$\omega_0 [t]$")
axs[1][0].set_ylabel(r"$(\omega_- - \omega_0) [t]$")
axs[2][0].set_ylabel(r"$W_0$")
legend = axs[0][1].legend(loc='upper center')
fig.tight_layout()

plt.savefig(f"plots/{os.path.basename(__file__).split('.')[0]}.pdf")
#plt.show()
