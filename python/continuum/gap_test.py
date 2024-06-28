import matplotlib.pyplot as plt
import numpy as np
import gzip

import __path_appender as __ap
__ap.append()
from extract_key import *
from create_zoom import *

MEV_FACTOR = 1e3

fig, ax = plt.subplots()

with gzip.open("data/continuum/test/gap.dat.gz", 'rt') as f_open:
    M = np.loadtxt(f_open)
M[0] = M[0] - extract_key("data/continuum/test/gap.dat.gz", "k_F")

for i in range(1, 4):
    M[i] *= MEV_FACTOR

ax.plot(M[0], (M[1] + M[2]), "k-", label=r"$\Delta_\mathrm{SC}$")
ax.plot(M[0], M[1], "--", label=r"$\Delta_\mathrm{Phonon}$")
ax.plot(M[0], M[2], "--", label=r"$\Delta_\mathrm{Coulomb}$")
ax.plot(M[0], M[3], "-", label=r"$\Delta_\mathrm{Fock}$")

create_zoom(ax, 0.1, 0.3, 0.3, 0.6, xlim=(-0.02, 0.02), ylim=(1.1 * np.min(M[2]), 1.05 * np.max(M[1])))

ax.set_xlabel(r"$k - k_\mathrm{F} [\sqrt{\mathrm{eV}}]$")
ax.set_ylabel(r"$\Delta [\mathrm{meV}]$")
ax.legend()
fig.tight_layout()

plt.show()