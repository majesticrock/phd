import matplotlib.pyplot as plt
import numpy as np
import gzip

plot_approx = False
folder = "test"

with gzip.open(f"data/continuum/{folder}/gap.dat.gz", 'rt') as f_open:
    M = np.loadtxt(f_open)
plt.plot(M[0], M[1], "-", label=r"exact $\theta$")

if plot_approx:
    with gzip.open(f"data/continuum/{folder}/gap_approx.dat.gz", 'rt') as f_open:
        M_approx = np.loadtxt(f_open)
    plt.plot(M_approx[0], M_approx[1], "-", label=r"approx. $\theta$")


plt.xlabel(r"$k [\sqrt{\mathrm{eV}}]$")
plt.ylabel(r"$\Delta [\mathrm{eV}]$")
plt.legend()
plt.tight_layout()
plt.show()