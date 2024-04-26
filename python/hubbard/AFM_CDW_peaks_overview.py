import numpy as np
import matplotlib.pyplot as plt
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib.continued_fraction as cf
from lib.iterate_containers import naming_scheme
from lib.extract_key import *
import lib.resolvent_peak as rp
# Calculates the resolvent in w^2

prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

Ts = np.array([0.])
Us = np.array([3.7, 3.725, 3.75, 3.8, 3.825, 3.85, 3.875, 3.9, 3.91,  # 10
                3.925, 3.94, 3.95, 3.96, 3.97, 3.975, 3.985, 3.99, 3.995, # 9
                3.997, 3.9985, 3.999, 3.9995, # 4
                4.0, # 1
                4.0005, 4.001, 4.0015, 4.003, # 4
                4.005, 4.01, 4.015, 4.025, 4.03, 4.04, 4.05, 4.06, 4.075, # 9
                4.09, 4.1, 4.125, 4.15, 4.175, 4.2, 4.25, 4.275, 4.3 # 10
                ])
Vs = np.array([1.])

folder = "data/modes/square/dos_3k/"
colors = ["orange", "purple"]

name_suffices = ["AFM", "CDW"]
fig, ax = plt.subplots()

for i, name_suffix in enumerate(name_suffices):
    peak_positions = np.zeros(len(Us))
    counter = 0
    for name in naming_scheme(Ts, Us, Vs):
        peak_positions[counter] = rp.Peak(f"{folder}{name}", name_suffix, (2, 3)).peak_position
        counter += 1

    u_data = (np.array([float(u) for u in Us]))#
    ax.plot(u_data, peak_positions, "X", label=f"{name_suffix}", color=colors[i])

ax.set_xlabel(r"$U / t$")
ax.set_ylabel(r"$z_0 / t$")
legend = plt.legend()

fig.tight_layout()


import os
plt.savefig(f"python/build/{os.path.basename(__file__).split('.')[0]}.pdf")
plt.show()
