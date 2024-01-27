import numpy as np
import matplotlib.pyplot as plt

fig, axs = plt.subplots(nrows=2, sharex=True, figsize=(6.4, 9))

xlims = (-5, 0)
ylims = (-0.5, 0)
alpha = .6
folder = "T0"

############################################################################
############                  Square lattice                    ############
############################################################################
import gzip
with gzip.open(f"../data/phases/square/{folder}/unkown_boundary.dat.gz", 'rt') as fp:
    u_boundary = np.loadtxt(fp)

axs[0].fill_between(u_boundary[0], u_boundary[1], np.full(len(u_boundary[0]), ylims[1]), label="$s$-wave", alpha=alpha)
axs[0].fill_between(u_boundary[0], u_boundary[1], np.full(len(u_boundary[0]), ylims[0]), label="New phases", alpha=alpha)
axs[0].fill_between([u_boundary[0][-1], xlims[1]], [ylims[1], ylims[1]], [ylims[0], ylims[0]], label=r"$\Delta$ too small", alpha=alpha)
axs[0].plot(u_boundary[0], u_boundary[1], "k-")
axs[0].axvline(u_boundary[0][-1], color="k", linestyle="--")

############################################################################
############                  Cubic lattice                     ############
############################################################################
with gzip.open(f"../data/phases/cube/{folder}/unkown_boundary.dat.gz", 'rt') as fp:
    u_boundary = np.loadtxt(fp)

axs[1].fill_between(u_boundary[0], u_boundary[1], np.full(len(u_boundary[0]), ylims[1]), label="$s$-wave", alpha=alpha)
axs[1].fill_between(u_boundary[0], u_boundary[1], np.full(len(u_boundary[0]), ylims[0]), label="New phases", alpha=alpha)
axs[1].fill_between([u_boundary[0][-1], xlims[1]], [ylims[1], ylims[1]], [ylims[0], ylims[0]], label=r"$\Delta$ too small", alpha=alpha)
axs[1].plot(u_boundary[0], u_boundary[1], "k-")
axs[1].axvline(u_boundary[0][-1], color="k", linestyle="--")

############################################################################
############                 Figure settings                    ############
############################################################################
axs[0].set_xlim(xlims[0], xlims[1])
axs[0].set_ylim(ylims[0], ylims[1])
axs[1].set_ylim(ylims[0], ylims[1])

axs[0].text(xlims[1] - np.abs(xlims[1] - xlims[0]) * 0.15, ylims[1] - np.abs(ylims[1] - ylims[0]) * 0.15, "(a)")
axs[1].text(xlims[1] - np.abs(xlims[1] - xlims[0]) * 0.15, ylims[1] - np.abs(ylims[1] - ylims[0]) * 0.15, "(b)")

axs[1].set_xlabel("$U / t$")
axs[0].set_ylabel("$V / t$")
axs[1].set_ylabel("$V / t$")

#axs[1].set_xticks([-2, -1, 0, 1, 2])
#axs[0].set_yticks([-2, -1, 0, 1, 2])
#axs[1].set_yticks([-2, -1, 0, 1, 2])

axs[0].legend(loc="upper left")
plt.tight_layout()
import os
plt.savefig(f"plots/{os.path.basename(__file__).split('.')[0]}.pdf")
#plt.show()
