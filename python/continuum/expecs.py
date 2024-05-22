import numpy as np
import matplotlib.pyplot as plt
import gzip

with gzip.open(f"data/continuum/test/expecs.dat.gz", 'rt') as f_open:
    M = np.loadtxt(f_open)
    
plt.plot(M[0], M[1], "-", label="$<n_k>$")
plt.plot(M[0], M[2], "-", label="$<f_k>$")

plt.legend()
plt.xlabel("$k$")
plt.ylabel("$<O>$")

plt.tight_layout()
plt.show()
