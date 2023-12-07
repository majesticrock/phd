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
Us = np.array([-2.0])
Vs = np.array([-0.1])

dos_based = ps.CURVEFAMILY(total_size(Ts, Us, Vs), True, ax)
momentum_based = ps.CURVEFAMILY(total_size(Ts, Us, Vs), True, ax)

momentum_based.set_shared_linestyle(":")
momentum_based.set_shared_kwargs(linewidth=2*plt.rcParams["lines.linewidth"])

momentum_based.set_individual_colors("nice2")
dos_based.set_individual_colors("nice2")

types = [["data/modes/square/dos_3k/", dos_based], ["data/modes/square/momentum_L=40/", momentum_based]]
name_suffix = "higgs_SC"

ax.set_yscale("log")

plot_upper_lim = 8.5

color_labels = []
for folder, curves in types:
    counter = 0
    for name in naming_scheme(Ts, Us, Vs):
        data_imag, data, w_lin, res = cf.resolvent_data(f"{folder}{name}", name_suffix, 0, plot_upper_lim, xp_basis=True)

        print(np.trapz(data_imag))
        if folder == types[0][0]:
            curves.plot(w_lin, 4 * data_imag)
        else:
            curves.plot(w_lin, data_imag)
        color_labels.append(f"$U={Us[counter]}$")
        counter += 1

#res.mark_continuum(ax)
linestyle_labels = ['DOS', 'Momentum']
color_and_linestyle_legends(ax, linestyle_labels=linestyle_labels, color_labels=color_labels, color_legend_title='', linestyle_legend_title="Method")

ax.set_xlabel(r"$\epsilon / t$")
ax.set_ylabel(r"$A(z)$")
fig.tight_layout()

import os
plt.savefig(f"python/build/{os.path.basename(__file__).split('.')[0]}.pdf")
plt.show()
