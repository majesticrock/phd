import numpy as np
import matplotlib.pyplot as plt
import gzip
import continued_fraction as cf
# Calculates the resolvent in w^2

prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

Ts = np.array([0.])
Us = np.array([-2.0])
Vs = np.array([-0.5])

folders = ["L=30", "L=40", "L=50", "L=60", "L=70"]
name_suffix = "SC"
fig, ax = plt.subplots()
types = [ "phase"]#"higgs",
lss = ["-", "-", "-", "--", ":", ]
plt.rcParams["lines.linewidth"] *= 2

ax.set_yscale("log")
axins = ax.inset_axes([0.1, 0.04, 0.2, 0.4])

plot_upper_lim = 8.5

for folder, ls in zip(folders, lss):
    for q, T in enumerate(Ts):
        for r, U in enumerate(Us):
            for s, V in enumerate(Vs):
                name = f"T={T}/U={U}_V={V}"
                w_vals = 20000
                w_lin = np.linspace(0, plot_upper_lim, w_vals, dtype=complex)
                w_lin += 1e-8j
                w_lin = w_lin**2

                for idx, type in enumerate(types):
                    res = cf.ContinuedFraction(f"data/{folder}/{name}", f"resolvent_{type}_{name_suffix}", True)

                    if(idx == 0):
                        ax.plot(np.sqrt(w_lin.real), -res.continued_fraction( np.copy(w_lin) ).imag, color=colors[q+r+s],
                            linestyle=ls, linewidth=(plt.rcParams["lines.linewidth"]+idx*2), label=f"$V={V}$")
                    else:
                        ax.plot(np.sqrt(w_lin.real), -res.continued_fraction( np.copy(w_lin) ).imag, color=colors[q+r+s],
                            linestyle=ls, linewidth=(plt.rcParams["lines.linewidth"]+idx*2))

res.mark_continuum(ax)
ax.legend()

ax.set_xlabel(r"$\epsilon / t$")
x1, x2, y1, y2 = 0, 0.15, 0.63, 650
axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)
axins.set_xticklabels([])
axins.set_yticklabels([])
ax.indicate_inset_zoom(axins, edgecolor="black")

fig.tight_layout()

import os
plt.savefig(f"python/build/{os.path.basename(__file__).split('.')[0]}.svg")
if(plt.rcParams["backend"] != "pgf"):
    plt.show()
