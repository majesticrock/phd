import numpy as np
import matplotlib.pyplot as plt
import lib.continued_fraction as cf
from lib.iterate_containers import *
import lib.plot_settings as ps
from lib.color_and_linestyle_legends import *
# Calculates the resolvent in w^2

prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']
fig, ax = plt.subplots()

Ts = np.array([0.])
Us = np.array([-2.0, -3.0])
Vs = np.array([-0.1])

dos_based = ps.CURVEFAMILY(total_size(Ts, Us, Vs), True, ax)
momentum_based = ps.CURVEFAMILY(total_size(Ts, Us, Vs), True, ax)

momentum_based.set_shared_linestyle(":")
momentum_based.set_shared_kwargs(linewidth=2*plt.rcParams["lines.linewidth"])

momentum_based.set_individual_colors("nice")
dos_based.set_individual_colors("nice")

types = [["data/modes/square/dos_disc=400/", dos_based], ["data/modes/square/momentum_L=20/", momentum_based]]
name_suffix = "higgs_sc"
element_names = ["a", "a+b", "a+ib"]

#ax.set_xscale("log")
ax.set_yscale("log")

plot_upper_lim = 8.5

color_labels = []
for folder, curves in types:
    for T, U, V in iterate_containers(Ts, Us, Vs):
        name = f"T={T}/U={U}_V={V}"

        w_vals = 20000
        w_lin = np.linspace(0, plot_upper_lim, w_vals, dtype=complex)
        w_lin += 1e-6j
        w_lin = w_lin**2
        data = np.zeros(w_vals)

        for idx, element in enumerate(element_names):
            file = f"{folder}{name}resolvent_{name_suffix}_{element}.dat.gz"
            res = cf.ContinuedFraction(f"{folder}{name}", f"resolvent_{name_suffix}_{element}", True)

            def dos(w):
                if idx==0:
                    return res.continued_fraction(w)
                else:
                    return np.sqrt(w) * res.continued_fraction(w)

            if idx == 0:
                data -= dos( np.copy(w_lin) ).imag
            elif idx == 1:
                data -= dos( np.copy(w_lin) ).imag
            elif idx == 2:
                data += dos( np.copy(w_lin) ).imag

        curves.plot(np.sqrt(w_lin.real), data)
        color_labels.append(f"$U={U}$")

#res.mark_continuum(ax)
linestyle_labels = ['DOS', 'Momentum']
color_and_linestyle_legends(ax, linestyle_labels=linestyle_labels, color_labels=color_labels, color_legend_title='', linestyle_legend_title="Method")

#legend = plt.legend()
#
#import matplotlib.lines as mlines
#dummy_lines = []
#dummy_lines.append(mlines.Line2D([],[], color="k", linestyle="-"))
#dummy_lines.append(mlines.Line2D([],[], color="k", linestyle="--", linewidth=2*plt.rcParams["lines.linewidth"]))
#legend_extra = plt.legend([dummy_lines[i] for i in [0,1]], [r"DOS", r"Momentum"], loc="upper center")
#
#ax.add_artist(legend)
#ax.add_artist(legend_extra)

ax.set_xlabel(r"$\epsilon / t$")
ax.set_ylabel(r"$A(z)$")
fig.tight_layout()

import os
plt.savefig(f"python/build/{os.path.basename(__file__).split('.')[0]}.pdf")
plt.show()
