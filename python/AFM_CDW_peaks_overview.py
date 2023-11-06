import numpy as np
import matplotlib.pyplot as plt
import lib.continued_fraction as cf
from lib.iterate_containers import iterate_containers
from lib.extract_key import *
# Calculates the resolvent in w^2

prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

Ts = np.array([0.])
Us = np.array([3.75, 3.8, 3.85, 3.9, 3.925, 3.95, 3.975, 3.99, 3.995, 3.999,
               4.001, 4.005, 4.01, 4.025, 4.05, 4.075, 4.1, 4.15, 4.2, 4.25])
Vs = np.array([1.])

folder = "data/modes/square/dos_900/"
element_names = ["a", "a+b", "a+ib"]
colors = ["orange", "purple"]

name_suffices = ["AFM", "CDW"]
fig, ax = plt.subplots()

for i, name_suffix in enumerate(name_suffices):
    peak_positions = np.zeros(len(Us))
    counter = 0
    for T, U, V in iterate_containers(Ts, Us, Vs):
        name = f"T={T}/U={U}/V={V}"
        data, data_real, w_lin, res = cf.resolvent_data(f"{folder}{name}", name_suffix, lower_edge=2, upper_edge=3, number_of_values=20000, imaginary_offset=1e-6, xp_basis=True, messages=False)

        peak_positions[counter] = w_lin[np.argmax(data)] #/ extract_key(f"{folder}{name}/resolvent_higgs_{name_suffix}.dat.gz", "Total Gap")
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
