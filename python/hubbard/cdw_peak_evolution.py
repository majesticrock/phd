import numpy as np
import matplotlib.pyplot as plt

import __path_appender as __ap
__ap.append()

import continued_fraction as cf
from iterate_containers import *
from extract_key import *
import resolvent_peak as rp
import plot_settings as ps



Ts = np.array([0.])
Us = np.array([-5.0, -4.5, -4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5])
Vs = np.array([1.2])

folders = ["./data/modes/square/dos_900/", "./data/modes/cube/dos_900/"]
element_names = ["a", "a+b", "a+ib"]

name_suffices = ["CDW"]
nrows=2
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
        
    weights = np.zeros(len(Us))
    peak_positions = np.zeros(len(Us))
    
    for name_suffix in name_suffices:
        counter = 0
        for T, U, V_a in iterate_containers(Ts, Us, Vs):
            V = 1.2 if i == 0 else 0.8
            name = f"T={T}/U={round(U, 4)}/V={V}"
            cont_edges = cf.continuum_edges(f"{folders[i]}{name}", f"higgs_{name_suffix}", True)
            lower = 0
            upper = cont_edges[0]

            peak_positions[counter], weights[counter] = rp.analyze_peak(f"{folders[i]}{name}", f"higgs_{name_suffix}", (lower, upper), peak_position_tol=1e-14, 
                                                                        reversed=True, **(peak_fit_params[1 if counter > len(Us) - 6 or counter < 5 else 0]))
            counter += 1
        
        label_subscript = name_suffix if name_suffix != "AFM" else "l.AFM"
        
        #us_same = np.array([np.concatenate((Us_square[::-1], Us_square)), np.concatenate((Us_cube[::-1], Us_cube))], dtype=object)
        plotters[0][i].plot(Us, peak_positions, label=f"$\\mathcal{{A}}_\\mathrm{{{label_subscript}}} (\\omega)$")
        plotters[1][i].plot(Us, np.exp(weights))

axs[0][0].set_title("Square lattice")
axs[0][1].set_title("Simple cubic lattice")

for i in range(2):
    axs[0][i].text(0.87, 0.5, f"(a.{i+1})", transform = axs[0][i].transAxes)
    axs[1][i].text(0.87, 0.5, f"(b.{i+1})", transform = axs[1][i].transAxes)

axs[nrows-1][0].set_xlabel(r"$U [t]$")
axs[nrows-1][1].set_xlabel(r"$U [t]$")
axs[0][0].set_ylabel(r"$\omega_0 [t]$")
axs[1][0].set_ylabel(r"$W_0$")
legend = axs[0][1].legend(loc='upper center')
fig.tight_layout()

plt.savefig(f"python/build/{os.path.basename(__file__).split('.')[0]}.pdf")
plt.show()
