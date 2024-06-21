import numpy as np
import matplotlib.pyplot as plt
import gzip
import __path_appender as __ap
__ap.append()

from color_and_linestyle_legends import *

fig, ax = plt.subplots()

with gzip.open(f"data/continuum/exact_theta/expecs.dat.gz", 'rt') as f_open:
    M = np.loadtxt(f_open)
    M[0] = M[0] - M[0][int(len(M[0]) / 2)]
    
ax.plot(M[0], M[1], ls="-", color="C0")
ax.plot(M[0], M[2], ls="-", color="C1")

with gzip.open(f"data/continuum/approx_theta/expecs.dat.gz", 'rt') as f_open:
    M = np.loadtxt(f_open)
    M[0] = M[0] - M[0][int(len(M[0]) / 2)]
    
ax.plot(M[0], M[1], ls="--", color="C0", linewidth=4)
ax.plot(M[0], M[2], ls="--", color="C1", linewidth=4)

with gzip.open(f"data/continuum/test/expecs.dat.gz", 'rt') as f_open:
    M = np.loadtxt(f_open)
    M[0] = M[0] - M[0][int(len(M[0]) / 2)]
    
ax.plot(M[0], M[1], ls=":", color="C0", linewidth=4)
ax.plot(M[0], M[2], ls=":", color="C1", linewidth=4)

color_and_linestyle_legends(ax, color_labels=[r"$\langle n_k \rangle$", r"$\langle f_k \rangle$"], 
                            linestyle_labels=["Exact interaction", "Approx. interaction", "With Coulomb"])

ax.set_xlabel(r"$k - k_\mathrm{F} [\sqrt{\mathrm{meV}}]$")
ax.set_ylabel("$<O>$")

fig.tight_layout()
plt.show()
