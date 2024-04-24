import matplotlib.pyplot as plt
import numpy as np
import gzip

with gzip.open("data/continuum/test/gap.dat.gz", 'rt') as f_open:
    M = np.loadtxt(f_open)
    
DISC = len(M)
def k_transform(u):
    return u / (1 - u)

k_space = k_transform(np.linspace(0., 1., DISC, endpoint=False))
plt.plot(k_space, M)
plt.xlabel(r"$k [\sqrt{meV}]$")
plt.ylabel(r"$\Delta [\mathrm{meV}]$")

plt.tight_layout()
plt.show()