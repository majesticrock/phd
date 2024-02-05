import numpy as np
import matplotlib.pyplot as plt

fig, axs = plt.subplots(nrows=2, sharex=True, figsize=(6.4, 9))

xlims = (-2, 2)
ylims = (-2, 2)
alpha = .6

############################################################################
############                  Square lattice                    ############
############################################################################
#Shade s-wave area
axs[0].fill_between([xlims[0], 0], [0, 0], [ylims[0], ylims[0]], alpha=alpha, label="$s$-wave")

#Shade CDW area
axs[0].fill_between(np.array([xlims[0], 0, xlims[1]]), np.array([0, 0, xlims[1] / 4]), np.array([ylims[1], ylims[1], ylims[1]]), alpha=alpha, label="CDW")

#Shade AFM area
axs[0].fill_between([0, xlims[1]], [0, xlims[1] / 4], [ylims[0], ylims[0]], alpha=alpha, label="AFM")

#Shade d-wave area
#axs[0].fill_between(micnas_d[0], micnas_d[1], -2 * np.ones(len(micnas_d[0])), color="red", alpha=.5, label="$d_{x^2 - y^2}$-wave")

#Plot boundaries
micnas_d = np.loadtxt("../data/micnas_d_wave.csv").transpose()
axs[0].plot(micnas_d[0], micnas_d[1], "k--", label="$d_{x^2 - y^2}$-boundary")
axs[0].plot(np.array([0, xlims[1]]), np.array([0, xlims[1]]) / 4,  "k-")
axs[0].plot(np.array([xlims[0], 0]), np.array([0, 0]),  "k-")
axs[0].plot(np.array([0, 0]),        np.array([0, ylims[0]]), "k-")

############################################################################
############                  Cubic lattice                     ############
############################################################################
#Shade s-wave area
axs[1].fill_between([xlims[0], 0], [0, 0], [ylims[0], ylims[0]], alpha=alpha, label="$s$-wave")

#Shade CDW area
axs[1].fill_between(np.array([xlims[0], 0, xlims[1]]), np.array([0, 0, xlims[1] / 6]), np.array([ylims[1], ylims[1], ylims[1]]), alpha=alpha, label="CDW")

#Shade AFM area
axs[1].fill_between([0, xlims[1]], [0, xlims[1] / 6], [ylims[0], ylims[0]], alpha=alpha, label="AFM")

#Plot boundaries
axs[1].plot(np.array([0, xlims[1]]), np.array([0, xlims[1]]) / 6,  "k-")
axs[1].plot(np.array([xlims[0], 0]), np.array([0, 0]),  "k-")
axs[1].plot(np.array([0, 0]),        np.array([0, ylims[0]]), "k-")

############################################################################
############                 Figure settings                    ############
############################################################################
axs[0].set_xlim(xlims[0], xlims[1])
axs[0].set_ylim(ylims[0], ylims[1])
axs[1].set_ylim(ylims[0], ylims[1])

axs[0].text(xlims[1] - np.abs(xlims[1] - xlims[0]) * 0.15, ylims[1] - np.abs(ylims[1] - ylims[0]) * 0.15, "(a)")
axs[1].text(xlims[1] - np.abs(xlims[1] - xlims[0]) * 0.15, ylims[1] - np.abs(ylims[1] - ylims[0]) * 0.15, "(b)")

axs[1].set_xlabel("$U [t]$")
axs[0].set_ylabel("$V [t]$")
axs[1].set_ylabel("$V [t]$")

#axs[1].set_xticks([-2, -1, 0, 1, 2])
#axs[0].set_yticks([-2, -1, 0, 1, 2])
#axs[1].set_yticks([-2, -1, 0, 1, 2])

axs[0].legend(loc="upper left")
plt.tight_layout()
import os
plt.savefig(f"plots/{os.path.basename(__file__).split('.')[0]}.pdf")
#plt.show()
