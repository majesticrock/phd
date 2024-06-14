import matplotlib.pyplot as plt
import numpy as np
import gzip

MEV_FACTOR = 1e3

fig, ax = plt.subplots()

with gzip.open(f"data/continuum/test/gap.dat.gz", 'rt') as f_open:
    M = np.loadtxt(f_open)

delta_abs = np.sqrt(M[1]**2 + M[2]**2)
delta_phase = np.arctan2(M[2], M[1])
ax.plot(M[0], MEV_FACTOR * delta_abs, "-", label=r"Abs SC")
ax.plot(M[0], MEV_FACTOR * M[3], "-", label=r"$\delta n$")

ax2 = ax.twinx()
ax2.plot(M[0], delta_phase / np.pi, "--", color="C3", label=r"Phase SC")
ax2.set_ylabel(r"arg[$\Delta$] / $\pi$")

ax.set_xlabel(r"$k [\sqrt{\mathrm{eV}}]$")
ax.set_ylabel(r"$\Delta [\mathrm{meV}]$")
fig.legend()
fig.tight_layout()

plt.show()