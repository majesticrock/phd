import matplotlib.pyplot as plt
import numpy as np

M = np.loadtxt("data/T0/U_modes/0.00.txt").transpose()
E = np.loadtxt("data/T0/U_modes/0.00_one_particle.txt").transpose()

for i in range(len(M)):
    plt.plot(M[i])

for i in range(len(E)):
    plt.plot(2*E[i], "--")

plt.xlabel(r"$k/\pi$")
plt.ylabel(r"$\epsilon / t$")
plt.show()