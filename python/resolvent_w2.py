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

folder = "data/modes/square/L=50/"
name_suffix = "SC"
fig, ax = plt.subplots()

#ax.set_xscale("log")
ax.set_yscale("log")

types = ["phase"]#"higgs", 
lss = ["-", "--", "-."]

plot_upper_lim = 8.5

for q, T in enumerate(Ts):
    for r, U in enumerate(Us):
        for s, V in enumerate(Vs):
            name = f"T={T}/U={U}_V={V}"
            w_vals = 20000
            w_lin = np.linspace(0, plot_upper_lim, w_vals, dtype=complex)
            w_lin += 1e-8j
            w_lin = w_lin**2
            
            for idx, type in enumerate(types):
                res = cf.ContinuedFraction(f"{folder}{name}", f"resolvent_{type}_{name_suffix}", True)

                if(idx == 0):
                    ax.plot(np.sqrt(w_lin.real), -res.continued_fraction( np.copy(w_lin) ).imag, color=colors[q+r+s],
                        linestyle=lss[idx], linewidth=(plt.rcParams["lines.linewidth"]+idx*2), label=f"$V={V}$")
                else:
                    ax.plot(np.sqrt(w_lin.real), -res.continued_fraction( np.copy(w_lin) ).imag, color=colors[q+r+s],
                        linestyle=lss[idx], linewidth=(plt.rcParams["lines.linewidth"]+idx*2))

res.mark_continuum(ax)
legend = plt.legend()

import matplotlib.lines as mlines
dummy_lines = []
dummy_lines.append(mlines.Line2D([],[], color="k", linestyle="-"))
dummy_lines.append(mlines.Line2D([],[], color="k", linestyle="--", linewidth=2*plt.rcParams["lines.linewidth"]))
legend_extra = plt.legend([dummy_lines[i] for i in [0,1]], [r"Amplitude", r"Phase"], loc="upper center")

ax.add_artist(legend)
ax.add_artist(legend_extra)

ax.set_xlim(-0.03, 8.5)
ax.set_xlabel(r"$\epsilon / t$")
ax.set_ylabel(r"$A(z)$")
fig.tight_layout()

import os
plt.savefig(f"python/build/{os.path.basename(__file__).split('.')[0]}.pdf")
plt.show()
