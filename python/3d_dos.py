import numpy as np
import matplotlib.pyplot as plt
import gzip

with gzip.open("data/3d_dos.dat.gz", 'rt') as f_open:
    half_dos = np.loadtxt(f_open)

dos = np.concatenate([half_dos[::-1], half_dos])
plt.plot(np.linspace(-3, 3, dos.size), dos)
plt.xlim(-3.1, 3.1)
plt.show()