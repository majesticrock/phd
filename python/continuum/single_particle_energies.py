import matplotlib.pyplot as plt
import numpy as np
import gzip

fig, ax = plt.subplots()

with gzip.open(f"data/continuum/test/one_particle_energies.dat.dz", 'rt') as f_open:
    M = np.loadtxt(f_open)

for i in range(len(M[1])):
    M[1][i] = -M[1][i]
    if(np.abs(M[1][i]) < 1e-5):
        break
ax.axhline(0, color="k")
ax.plot(M[0], M[1], "-", label="With Coulomb")

def bare(k):
    return 0.5 * k * k - 9.3
ax.plot(M[0], bare(M[0]), "--", label="Bare")


ax.set_xlabel(r"$k [\sqrt{\mathrm{eV}}]$")
ax.set_ylabel(r"$E(k) [\mathrm{eV}]$")
ax.legend()
fig.tight_layout()

plt.show()