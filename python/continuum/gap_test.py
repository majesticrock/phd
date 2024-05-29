import matplotlib.pyplot as plt
import numpy as np
import gzip

MEV_FACTOR = 1e3

fig, ax = plt.subplots()

with gzip.open(f"data/continuum/test/gap.dat.gz", 'rt') as f_open:
    M = np.loadtxt(f_open)
ax.plot(M[0], MEV_FACTOR * M[1], "-", label=r"SC")
ax.plot(M[0], MEV_FACTOR * M[2], "-", label=r"OCC")

ax.set_xlabel(r"$k [\sqrt{\mathrm{eV}}]$")
ax.set_ylabel(r"$\Delta [\mathrm{meV}]$")
ax.legend()
fig.tight_layout()

plt.show()