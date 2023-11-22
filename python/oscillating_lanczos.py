import numpy as np
import matplotlib.pyplot as plt
import gzip

T = 0.
U = -2.0
V = -0.1

use_XP = True

folder = "data/modes/cube/dos_6k/"
name_suffix = "phase_SC"
name = f"T={T}/U={U}/V={V}/"
element_names = ["a", "a+b", "a+ib"]

file = f"{folder}{name}one_particle.dat.gz"
with gzip.open(file, 'rt') as f_open:
    one_particle = np.abs(np.loadtxt(f_open).flatten())
    roots = np.array([np.min(one_particle) * 2, np.max(one_particle) * 2])**2
    a_inf = (roots[0] + roots[1]) * 0.5
    b_inf = ((roots[1] - roots[0]) * 0.25)

if use_XP:
    file = f"{folder}{name}resolvent_{name_suffix}.dat.gz"
else:
    element = element_names[0]
    file = f"{folder}{name}resolvent_{name_suffix}_{element}.dat.gz"
with gzip.open(file, 'rt') as f_open:
    M = np.loadtxt(f_open)
    A = M[0]
    B = M[1]

fig, ax = plt.subplots()

osc_frequency = 6
A_star = np.zeros(len(A))
for i in range(len(A)):
    if len(A) - i < osc_frequency:
        A_star[i] = 0
    else:
        for j in range(osc_frequency):
            A_star[i] += A[i + j]
        A_star[i] /= osc_frequency
        A_star[i] -= A[i]
    
A_star = A_star[1:-osc_frequency]
n_space = np.linspace(1, len(A_star), len(A_star))
ax.plot(n_space, A_star, ls="-", marker='x', label="$a_i^*$")

from scipy.optimize import curve_fit
def func(n, a, alpha, phi0):
    return a * n**(-1-alpha) * np.sin((2*n - 1) + 2*phi0)

fit_cut = 3
popt, pcov = curve_fit(func, n_space[fit_cut:], A_star[fit_cut:])
ax.plot(np.linspace(1, len(A_star), 500), func(np.linspace(1, len(A_star), 500), *popt))
print(popt)

ax.legend()
ax.set_xlabel("Iteration $i$")
ax.set_ylabel("Lanczos coefficient")
#ax.set_ylim(13, 34)
fig.tight_layout()

import os
plt.savefig(f"python/build/{os.path.basename(__file__).split('.')[0]}.pdf")
plt.show()
