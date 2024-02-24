import numpy as np
import matplotlib.pyplot as plt
import gzip

T = 0.
U = -2.5
V = 0.0

use_XP = True

folder = "data/modes/square/test/"
name_suffix = "higgs_SC"
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

#from scipy.signal import find_peaks
#first_artifact = find_peaks(B[1:], prominence=5e2, width=1)
#print(first_artifact)

fig, ax = plt.subplots()
ax.plot(A, ls="-", marker='x', label="$a_i$")
ax.plot(np.sqrt(B), ls="-", marker='o', label="$b_i$")
ax.axhline(a_inf, linestyle="-" , color="k", label="$a_\\infty$")
ax.axhline(b_inf, linestyle="--", color="k", label="$b_\\infty$")
ax.legend()
ax.set_xlabel("Iteration $i$")
ax.set_ylabel("Lanczos coefficient")
#ax.set_ylim(13, 34)
fig.tight_layout()

import os
plt.savefig(f"python/build/{os.path.basename(__file__).split('.')[0]}.pdf")
plt.show()
