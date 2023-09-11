import numpy as np
import matplotlib.pyplot as plt
import gzip

file = "data/modes/square/test/T=0.0/U=-2.0_V=-0.1/one_particle.dat.gz"
with gzip.open(file, 'rt') as f_open:
    M = np.loadtxt(f_open).transpose()

plt.plot(M)

plt.ylabel(r"$\epsilon / t$")
plt.xlabel(r"$k / \pi$")
plt.tight_layout()

import os
plt.savefig(f"python/build/{os.path.basename(__file__).split('.')[0]}.pdf")
plt.show()