import numpy as np
import matplotlib.pyplot as plt

fig, axs = plt.subplots(nrows=2, sharex=True, figsize=(6.4, 9))

############################################################################
############                  Square lattice                    ############
############################################################################
#Shade s-wave area
axs[0].fill_between([-2, 0], [0, 0], [-2, -2], color="blue", alpha=.5, label="$s$-wave")

#Shade AFM area
axs[0].fill_between([0, 2], [0, .5], [-2, -2], color="green", alpha=.5, label="AFM")

#Shade CDW area
axs[0].fill_between(np.array([-2, 0, 2]), np.array([0, 0, .5]), np.array([2, 2, 2]), color="orange", alpha=.5, label="CDW")

#Shade d-wave area
#axs[0].fill_between(micnas_d[0], micnas_d[1], -2 * np.ones(len(micnas_d[0])), color="red", alpha=.5, label="$d_{x^2 - y^2}$-wave")

#Plot boundaries
micnas_d = np.loadtxt("../data/micnas_d_wave.csv").transpose()
axs[0].plot(micnas_d[0], micnas_d[1], "k--", label="$d_{x^2 - y^2}$-boundary")
axs[0].plot(np.array([0, 2]), 0.25 * np.array([0, 2]),  "k-")
axs[0].plot(np.array([-2, 0]),       np.array([0, 0]),  "k-")
axs[0].plot(np.array([0, 0]),        np.array([0, -2]), "k-")

############################################################################
############                  Cubic lattice                     ############
############################################################################
#Shade s-wave area
axs[1].fill_between([-2, 0], [0, 0], [-2, -2], color="blue", alpha=.5, label="$s$-wave")

#Shade AFM area
axs[1].fill_between([0, 2], [0, 1./3.], [-2, -2], color="green", alpha=.5, label="AFM")

#Shade CDW area
axs[1].fill_between(np.array([-2, 0, 2]), np.array([0, 0, 1./3.]), np.array([2, 2, 2]), color="orange", alpha=.5, label="CDW")

#Plot boundaries
axs[1].plot(np.array([0, 2]),  np.array([0, 2]) / 6., "k-")
axs[1].plot(np.array([-2, 0]), np.array([0, 0]),      "k-")
axs[1].plot(np.array([0, 0]),  np.array([0, -2]),     "k-")

############################################################################
############                 Figure settings                    ############
############################################################################
axs[0].set_xlim(-2, 2)
axs[0].set_ylim(-2, 2)
axs[1].set_ylim(-2, 2)

axs[0].text(1.7, 1.7, "(a)")
axs[1].text(1.7, 1.7, "(b)")

axs[1].set_xlabel("$U / t$")
axs[0].set_ylabel("$V / t$")
axs[1].set_ylabel("$V / t$")

axs[1].set_xticks([-2, -1, 0, 1, 2])
axs[0].set_yticks([-2, -1, 0, 1, 2])
axs[1].set_yticks([-2, -1, 0, 1, 2])

axs[0].legend(loc="upper left")
plt.tight_layout()
plt.savefig("plots/phase_diagram.pdf")