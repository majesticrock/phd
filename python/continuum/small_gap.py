import numpy as np
import matplotlib.pyplot as plt
import gzip


def analytical(U, w_D, eps_F):
    rho_0 = 3 / (4 * eps_F)
    return 2 * w_D * np.exp(- 1 / (rho_0 * U))

folder = "test"
with gzip.open(f"data/continuum/{folder}/small_U_gap.dat.gz", 'rt') as f_open:
    M = np.loadtxt(f_open).transpose()
    
plt.plot(M[0], 1e3 * M[1], label="Self-consistency")
#popt, pcov = curve_fit(analytical, M[0], M[1], p0=(30., 3. / 4000.))

plt.plot(M[0], 1e3 * analytical(M[0], 0.01, 9.3), "--", label="Analytical prediction")

plt.legend()
plt.xlabel(r"$g$ [eV]")
plt.ylabel(r"$\Delta$ [meV]")

plt.tight_layout()
plt.show()