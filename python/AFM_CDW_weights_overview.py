import numpy as np
import matplotlib.pyplot as plt
from lib.iterate_containers import naming_scheme
from lib.extract_key import *
import lib.resolvent_peak as rp

prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

Ts = np.array([0.])
Us = np.array([3.7, 3.725, 3.75, 3.775, 3.8, 3.825, 3.85, 3.875, 3.9, 3.91,  # 10
               3.925, 3.94, 3.95, 3.96, 3.97, 3.975, 3.985, 3.99, 3.995, # 9
               3.997, 3.9985, 3.999, 3.9995, # 4
               4.0, # 1
               4.0005, 4.001, 4.0015, 4.003, # 4
               4.005, 4.01, 4.015, 4.025, 4.03, 4.04, 4.05, 4.06, 4.075, # 9
               4.09, 4.1, 4.125, 4.15, 4.175, 4.2, 4.225, 4.25, 4.275, 4.3 # 10
               ])
Vs = np.array([1.])

folder = "data/modes/square/dos_900/"
name_suffices = ["AFM", "CDW"]
fig, ax = plt.subplots()
u_data = np.array([float(u) for u in Us])

for i, name_suffix in enumerate(name_suffices):
    weights = np.zeros(len(Us))
    cont_weights = np.zeros(len(Us))
    counter = 0
    for name in naming_scheme(Ts, Us, Vs):
        peak = rp.Peak(f"{folder}{name}", name_suffix, initial_search_bounds=(1., 4.))
        peak.improved_peak_position(x0_offset=1e-2, gradient_epsilon=1e-12)
        popt, pcov, w_space, y_data = peak.fit_real_part(range=0.01, begin_offset=1e-10, reversed=True)

        weights[counter] = np.pi * np.exp(popt[1])
        cont_weights[counter] = peak.resolvent.weight_of_continuum(2000, 0)
        counter += 1

    ax.plot(u_data, weights, marker="X", ls="-", label=f"{name_suffix}")
    ax.plot(u_data, cont_weights, marker="v", ls="-", label=f"{name_suffix} - Cont")

ax.set_xlabel(r"$U / t$")
ax.set_ylabel(r"$w_0 \cdot t$")
ax.legend()
fig.tight_layout()

import os
plt.savefig(f"python/build/{os.path.basename(__file__).split('.')[0]}.pdf")
plt.show()
