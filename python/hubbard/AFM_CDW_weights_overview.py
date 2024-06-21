import numpy as np
import matplotlib.pyplot as plt
import __path_appender as __ap
__ap.append()
from iterate_containers import naming_scheme
from extract_key import *
import resolvent_peak as rp
import continued_fraction as cf

prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

Ts = np.array([0.])
Us = np.array([ -0.3, -0.275, -0.25, -0.2, -0.175, -0.15, -0.125, -0.1, -0.09, -0.075, -0.06, -0.05, -0.04, -0.03, -0.025, -0.01, -0.01, -0.005, -0.003, -0.0015, -0.001, -0.0005, 
                0., 
                0.0005, 0.00, 0.0015, 0.003, 0.005, 0.01, 0.015, 0.025, 0.03, 0.04, 0.05, 0.06, 0.075, 0.09, 0.1, 0.125, 0.15, 0.175, 0.2, 0.25, 0.275, 0.3
            ]) 
Vs = np.array([1.])

square = False
Us += (4 if square else 6)
folder = "data/modes/" + ("square" if square else "cube") + "/dos_3k/"
name_suffices = ["AFM", "CDW"]
fig, ax = plt.subplots()
u_data = np.array([float(u) for u in Us])

for i, name_suffix in enumerate(name_suffices):
    weights = np.zeros(len(Us))
    cont_weights = np.zeros(len(Us))
    counter = 0
    for name in naming_scheme(Ts, Us, Vs):
        peak = rp.Peak(f"{folder}{name}", name_suffix, initial_search_bounds=(1., cf.continuum_edges(f"{folder}{name}", name_suffix, xp_basis=True)[0]))
        peak.improved_peak_position(xtol=1e-13)
        popt, pcov, w_space, y_data = peak.fit_real_part(range=0.01, begin_offset=1e-10, reversed=True)

        weights[counter] = np.exp(popt[1])
        #cont_weights[counter] = peak.resolvent.weight_of_continuum(2000, 0)
        counter += 1

    ax.plot(u_data, weights, marker="X", ls="-", label=f"{name_suffix}")
    #ax.plot(u_data, cont_weights, marker="v", ls="-", label=f"{name_suffix} - Cont")

ax.set_xlabel(r"$U / t$")
ax.set_ylabel(r"$w_0 \cdot t$")
ax.legend()
fig.tight_layout()

import os
plt.savefig(f"python/build/{os.path.basename(__file__).split('.')[0]}.pdf")
plt.show()
