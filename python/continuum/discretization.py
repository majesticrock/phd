import matplotlib.pyplot as plt
import numpy as np
import gzip

MEV_FACTOR = 1e3

fig, ax = plt.subplots()

with gzip.open(f"data/continuum/test/gap.dat.gz", 'rt') as f_open:
    M = np.loadtxt(f_open)
M[0] = M[0] / M[0][int(0.5 * len(M[0])) - 1]

ax.plot(M[0], label="Discretization")

ax.set_xlabel(r"$n$")
ax.set_ylabel(r"$k / k_\mathrm{F}$")
fig.legend()
fig.tight_layout()

plt.show()