import numpy as np
import matplotlib.pyplot as plt

fig, axs = plt.subplots(nrows=2, sharex=True, figsize=(6.4, 9))
folder = "u_bound_1"

xlims = (-4, 4)
ylims = (-2, 2)
alpha = .6

############################################################################
############                  Square lattice                    ############
############################################################################
import gzip
with gzip.open(f"../data/phases/square/{folder}/unkown_boundary.dat.gz", 'rt') as fp:
    u_boundary = np.loadtxt(fp)

#Shade s-wave area
axs[0].fill_between(u_boundary[0], u_boundary[1], np.full(len(u_boundary[0]), 0), label="$s$-wave", alpha=alpha)

#Shade CDW area
axs[0].fill_between(np.array([xlims[0], 0, xlims[1]]), np.array([0, 0, xlims[1] / 4]), np.array([ylims[1], ylims[1], ylims[1]]), alpha=alpha, label="CDW")

#Shade AFM area
axs[0].fill_between([0, xlims[1]], [0, xlims[1] / 4], [ylims[0], ylims[0]], alpha=alpha, label="AFM")

# unknown regions
axs[0].fill_between(u_boundary[0], u_boundary[1], np.full(len(u_boundary[0]), ylims[0]), label="Likely PS", alpha=alpha)
axs[0].fill_between([u_boundary[0][-1], 0], [0, 0], [ylims[0], ylims[0]], label=r"$\Delta$ too small", 
                    alpha=alpha, facecolor="C0", edgecolor="black", hatch="///", linewidth=2.0)
axs[0].plot(u_boundary[0], u_boundary[1], "k-")
#axs[0].axvline(u_boundary[0][-1], ymin=0, ymax=0.5, color="k", linestyle=":")

#Plot boundaries
micnas_d = np.loadtxt("../data/micnas_d_wave.csv").transpose()
axs[0].plot(micnas_d[0], micnas_d[1], "k--", label="$d_{x^2 - y^2}$-boundary")
axs[0].plot(np.array([0, xlims[1]]), np.array([0, xlims[1]]) / 4,  "k-")
axs[0].plot(np.array([xlims[0], 0]), np.array([0, 0]),  "k-")
axs[0].plot(np.array([0, 0]),        np.array([0, ylims[0]]), "k-")

############################################################################
############                  Cubic lattice                     ############
############################################################################
with gzip.open(f"../data/phases/cube/{folder}/unkown_boundary.dat.gz", 'rt') as fp:
    u_boundary = np.loadtxt(fp)

#Shade s-wave area
axs[1].fill_between(u_boundary[0], u_boundary[1], np.full(len(u_boundary[0]), 0), label="$s$-wave", alpha=alpha)

#Shade CDW area
axs[1].fill_between(np.array([xlims[0], 0, xlims[1]]), np.array([0, 0, xlims[1] / 6]), np.array([ylims[1], ylims[1], ylims[1]]), alpha=alpha, label="CDW")

#Shade AFM area
axs[1].fill_between([0, xlims[1]], [0, xlims[1] / 6], [ylims[0], ylims[0]], alpha=alpha, label="AFM")

# unknown regions
axs[1].fill_between(u_boundary[0], u_boundary[1], np.full(len(u_boundary[0]), ylims[0]), label="Likely PS", alpha=alpha)
axs[1].fill_between([u_boundary[0][-1], 0], [0, 0], [ylims[0], ylims[0]], label=r"$\Delta$ too small", 
                    alpha=alpha, facecolor="C0", edgecolor="black", hatch="///", linewidth=2.0)
axs[1].plot(u_boundary[0], u_boundary[1], "k-")
#axs[1].axvline(u_boundary[0][-1], ymin=0, ymax=0.5, color="k", linestyle="-")

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

axs[0].text(-3.85, 1.65, "(a)")
axs[1].text(-3.85, 1.65, "(b)")

axs[1].set_xlabel("$U [t]$")
axs[0].set_ylabel("$V [t]$")
axs[1].set_ylabel("$V [t]$")

axs[0].legend(loc="upper right", ncol=2, columnspacing=1)
plt.tight_layout()
import os
plt.savefig(f"plots/{os.path.basename(__file__).split('.')[0]}.pdf")
#plt.show()
