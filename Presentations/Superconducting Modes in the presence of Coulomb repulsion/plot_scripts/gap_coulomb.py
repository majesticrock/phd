import matplotlib.pyplot as plt
import __path_appender as __ap
__ap.append()
from get_data import *
import os
import string

fig, axes = plt.subplots(nrows=2, ncols=2, sharex="col", figsize=(6.4, 6.4))

for i, (axis, screening) in enumerate(zip(axes, [1e-4, 1])):
    main_df = load_panda("continuum", "offset_20", "gap.json.gz", 
                        **continuum_params(N_k=20000, T=0, coulomb_scaling=1, screening=screening, k_F=4.25, g=0.8, omega_D=10))
    pd_data = main_df["data"]
    pd_data["ks"] /= main_df["k_F"]

    if pd_data["Delta_Coulomb"][0] > 0:
        pd_data["Delta_Phonon"] *= -1
        pd_data["Delta_Coulomb"] *= -1
    
    for j, ax in enumerate(axis):
        ax.plot(pd_data["ks"], pd_data["Delta_Phonon"] + pd_data["Delta_Coulomb"], "k-", label=r"$\Delta$")
        pd_data.plot(x="ks", y=["Delta_Phonon", "Delta_Coulomb", "Delta_Fock"], ax=ax, style=['--', '--', ':'], 
                    label=[r"$\Delta_\mathrm{Ph}$", r"$\Delta_\mathrm{C}$", r"$\epsilon_\mathrm{C}$"], legend=False)
        ax.text(0.015, 0.88, f"({string.ascii_lowercase[i]}.{j+1})", transform=ax.transAxes)
    
    axis[0].set_ylabel(r"$\Delta [\mathrm{meV}]$")
    axis[1].yaxis.tick_right()
    axis[1].yaxis.set_ticks_position('both')
    axis[0].set_xlim(1-0.004, 1.004)

axes[0][1].set_ylim(-0.45, 0.02)
axes[1][1].set_ylim(-0.45, 0.02)
axes[1][0].set_xlabel(r"$k / k_\mathrm{F}$")
axes[1][1].set_xlabel(r"$k / k_\mathrm{F}$")
axes[1][1].set_xticks([0, 1, 2])

legend = axes[0][1].legend(handlelength=1.25, loc = 'lower right')
bb = legend.get_bbox_to_anchor().transformed(axes[0][1].transAxes.inverted())
SHIFT = 0.025
bb.x0 += SHIFT
bb.x1 += SHIFT
bb.y0 -= SHIFT
bb.y1 -= SHIFT
legend.set_bbox_to_anchor(bb, transform=axes[0][1].transAxes)

fig.subplots_adjust(hspace=0.05, wspace=0.05)
fig.savefig(f"plots/{os.path.basename(__file__).split('.')[0]}.pdf")