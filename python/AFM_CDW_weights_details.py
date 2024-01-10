import numpy as np
import matplotlib.pyplot as plt
from lib.iterate_containers import naming_scheme
from lib.extract_key import *
import lib.resolvent_peak as rp
import lib.continued_fraction as cf

prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

Ts = np.array([0.])
Us = np.array([ -0.3, -0.2, -0.15, -0.1, -0.07, -0.05, -0.04, -0.03, -0.02, -0.015, -0.01, 
                -0.007, -0.005, -0.004, -0.003, -0.002, -0.0015, -0.001, -0.0007, 
                -0.0005, -0.0004, -0.0003, -0.0002, -0.00015, -0.0001, 
                0.0001, 0.00015, 0.0002, 0.0003, 0.0004, 0.0005, 
                0.0007, 0.001, 0.0015, 0.002, 0.003, 0.004, 0.005, 0.007, 
                0.01, 0.015, 0.02, 0.03, 0.04, 0.05, 0.07, 0.1, 0.15, 0.2, 0.3]) 
Vs = np.array([1.])
u_data = (np.array([(float(u)) for u in Us]))

square = True
Us += (4 if square else 6)
folder = "data/modes/" + ("square" if square else "cube") + "/dos_3k/"

name_suffix = "CDW"
fig, ax = plt.subplots()


weights = np.zeros(len(Us))
counter = 0
    
for name in naming_scheme(Ts, Us, Vs):
    peak = rp.Peak(f"{folder}{name}", name_suffix, initial_search_bounds=(1., cf.continuum_edges(f"{folder}{name}", name_suffix, xp_basis=True)[0]))
    
    peak.improved_peak_position(xtol=1e-13)
    popt, pcov, w_space, y_data = peak.fit_real_part(range=0.01, begin_offset=1e-10, reversed=True)

    weights[counter] = (popt[1])
    counter += 1

ax.plot(u_data, weights, marker="X", ls="-", label=f"{name_suffix}")

ax.set_xlabel(r"$U / t$")
ax.set_ylabel(r"$\ln(w_0 \cdot t)$")
ax.legend()
fig.tight_layout()

import os
plt.savefig(f"python/build/{os.path.basename(__file__).split('.')[0]}.pdf")
plt.show()
