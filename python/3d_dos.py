import numpy as np
import matplotlib.pyplot as plt
import gzip

with gzip.open("data/2d_dos.dat.gz", 'rt') as f_open:
    half_dos = np.loadtxt(f_open)

num_zero = int(half_dos.size * 0.55)
dos = np.concatenate([np.zeros(num_zero), half_dos[::-1], half_dos, np.zeros(num_zero)])
dos /= np.max(dos)
plt.plot(np.linspace(-3.1, 3.1, dos.size), dos, label="2D")


with gzip.open("data/3d_dos.dat.gz", 'rt') as f_open:
    half_dos = np.loadtxt(f_open)
    
num_zero = int(half_dos.size * 0.05 * (2. / 3.))
dos = np.concatenate([np.zeros(num_zero), half_dos[::-1], half_dos, np.zeros(num_zero)])
dos /= np.max(dos)
plt.plot(np.linspace(-3.1, 3.1, dos.size), dos, label="3D")

#from scipy.special import *
#off = 1e-2
#x = np.linspace(-1+off, 1-off, dos.size)
#def func(x, y):
#    return ellipk(1 - ((x+y)**2)/4) / np.sqrt(1 - x*x)
#
#plt.plot(x, func(x, 1))

plt.xlabel("$\\gamma$ / $t$")
plt.ylabel("$\\rho(\\gamma)$ / a.u.")
plt.legend()
plt.tight_layout()
plt.show()