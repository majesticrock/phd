import numpy as np
import matplotlib.pyplot as plt
import gzip

with gzip.open("data/2d_dos.dat.gz", 'rt') as f_open:
    M = np.loadtxt(f_open)
half_dos = M[1]
half_abscissa = M[0]

dos = np.concatenate([[0.0, 0.0], half_dos[::-1], half_dos, [0.0, 0.0]])
dos /= np.max(dos)
abscissa = np.concatenate([[-3.1, -half_abscissa[-1]], -half_abscissa[::-1], half_abscissa, [half_abscissa[-1], 3.1]])
plt.plot(abscissa, dos, label="2D")


with gzip.open("data/3d_dos.dat.gz", 'rt') as f_open:
    M = np.loadtxt(f_open)
half_dos = M[1]
half_abscissa = M[0]
    
dos = np.concatenate([[0.0, 0.0], half_dos[::-1], half_dos, [0.0, 0.0]])
dos /= np.max(dos)
abscissa = np.concatenate([[-3.1, -half_abscissa[-1]], -half_abscissa[::-1], half_abscissa, [half_abscissa[-1], 3.1]])
plt.plot(abscissa, dos, label="3D")

plt.xlabel("$\\gamma$ / $t$")
plt.ylabel("$\\rho(\\gamma)$ / a.u.")
plt.legend()
plt.tight_layout()
plt.show()