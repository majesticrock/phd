import numpy as np
import matplotlib.pyplot as plt

M = np.loadtxt("data/deviations_hermiticity.csv", delimiter=";").transpose()

plt.plot(M[0], M[1], "-o")
plt.yscale("log")

plt.text(10, 0.1, r"$\epsilon_0(k) = -\frac{U}{2L^2} \cdot \sum_q \frac{\langle g_q \rangle}{\langle g_k \rangle} (1 - 2\langle n_k \rangle)$", fontsize=20)
plt.xlabel(r"$L$")
plt.ylabel(r"Deviation")
plt.tight_layout()

import os
plt.savefig(f"python/build/{os.path.basename(__file__).split('.')[0]}.pdf")
plt.show()
