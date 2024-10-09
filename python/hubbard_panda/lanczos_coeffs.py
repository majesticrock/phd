import numpy as np
import matplotlib.pyplot as plt
import __path_appender as __ap
__ap.append()

import continuum_boundaries as cb

T = 0.
U = -2.5
V = -0.1

from get_data import load_panda, hubbard_params
pd_data = load_panda("hubbard/square", "test", "dispersions.json.gz", **hubbard_params(0.0, -2.5, 0.0))

kx = 20
ky = 2
N_total = pd_data["discretization"]
n_k = 22

single = cb.SingleParticle(pd_data["gap_parameters"])
lower_bound, upper_bound = cb.continuum_bounds(single.dispersion, kx / N_total, ky / N_total, N_total)

print(lower_bound, upper_bound)

a_inf = (upper_bound**2 + lower_bound**2) * 0.5
b_inf = (upper_bound**2 - lower_bound**2) * 0.25

A = pd_data["resolvents.phase_SC_a"][n_k]["a_i"]
B = pd_data["resolvents.phase_SC_a"][n_k]["b_i"]

fig, ax = plt.subplots()
ax.plot(A, ls="-", marker='x', label="$a_i$")
ax.plot(np.sqrt(B), ls="-", marker='o', label="$b_i$")
ax.axhline(a_inf, linestyle="-" , color="k", label="$a_\\infty$")
ax.axhline(b_inf, linestyle="--", color="k", label="$b_\\infty$")
ax.legend()
ax.set_xlabel("Iteration $i$")
ax.set_ylabel("Lanczos coefficient")
#ax.set_ylim(13, 34)
fig.tight_layout()

import os
plt.savefig(f"python/build/{os.path.basename(__file__).split('.')[0]}.pdf")
plt.show()
