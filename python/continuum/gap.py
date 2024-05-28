import matplotlib.pyplot as plt
import numpy as np
import gzip

MEV_FACTOR = 1e3

def energy(k, delta):
    return np.sqrt((0.5 * k**2 - 9.3)**2 + np.abs(delta)**2)

fig, ax = plt.subplots()

with gzip.open(f"data/continuum/exact_theta/gap.dat.gz", 'rt') as f_open:
    M = np.loadtxt(f_open)
ax.plot(M[0], MEV_FACTOR * M[1], "-", label=r"$\theta(\omega_\mathrm{D} - |\epsilon_k - \epsilon_l|)$")

ax2 = ax.twinx()
ax2.plot(M[0], energy(M[0], M[1]), "--")

with gzip.open(f"data/continuum/approx_theta/gap.dat.gz", 'rt') as f_open:
    M_approx = np.loadtxt(f_open)
ax.plot(M_approx[0], MEV_FACTOR * M_approx[1], "-", label=r"$\theta(\omega_\mathrm{D} - \epsilon_k) \theta(\omega_\mathrm{D} - \epsilon_l)$")
ax2.plot(M_approx[0], energy(M_approx[0], M_approx[1]), "--")

ax.set_xlabel(r"$k [\sqrt{\mathrm{eV}}]$")
ax.set_ylabel(r"$\Delta [\mathrm{meV}]$")
ax2.set_ylabel(r"$E(k) [\mathrm{eV}]$")
ax.legend()
fig.tight_layout()

plt.show()