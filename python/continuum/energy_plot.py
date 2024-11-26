import matplotlib.pyplot as plt
import numpy as np
import __path_appender as __ap
__ap.append()
from get_data import *
import os

fig, ax = plt.subplots()
def energy(xi, delta):
    return np.sqrt(xi**2 + 1e-6 * delta**2)

GS = [1.15, 1.175, 1.2]
SCREENING = 0
OFFS = 0.0001

for g in GS:
    main_df = load_panda("continuum", "offset_20", "gap.json.gz", print_date=False,
                        **continuum_params(N_k=20000, T=0, coulomb_scaling=int(SCREENING!=0), screening=SCREENING, k_F=4.25, g=g, omega_D=10))
    pd_data = main_df["data"]
    pd_data["ks"] /= main_df["k_F"]
    if pd_data["Delta_Coulomb"][0] > 0:
        pd_data["Delta_Phonon"] *= -1
        pd_data["Delta_Coulomb"] *= -1
    plot_data = pd_data.query(f"ks > {1-OFFS} & ks < {1+OFFS}")
    ax.plot(plot_data["ks"], 
            1e3 * energy(plot_data["xis"] + 1e-3 * plot_data["Delta_Fock"],
                         plot_data["Delta_Phonon"] + plot_data["Delta_Coulomb"]) / main_df["Delta_max"] - 1, 
            label=f"$g={g}$")
    ax.set_xlabel(r"$k / k_\mathrm{F}$")

ax.set_ylabel(r"$E / \Delta_\mathrm{max} - 1$")
ax.legend(loc="upper right", ncols=2)
#ax.set_ylim(-1, 4)
fig.subplots_adjust(wspace=0.2, hspace=0.1)
plt.show()