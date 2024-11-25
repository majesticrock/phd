import matplotlib.pyplot as plt
import numpy as np
import __path_appender as __ap
__ap.append()
from get_data import *
import os

fig, axes = plt.subplots(ncols=3, sharey=True)
def energy(xi, delta):
    return np.sqrt(xi**2 + 1e-6 * delta**2)

GS = [2.15, 2.175, 2.2]
SCREENINGS = [0, 1, 1]
OFFS = [0.0014, 0.001, 0.0005]

for i in range(3):
    for g in GS:
        main_df = load_panda("continuum", "offset_25", "gap.json.gz", print_date=False,
                            **continuum_params(N_k=30000, T=0, coulomb_scaling=float(i != 0), screening=SCREENINGS[i], k_F=4.25, g=g, omega_D=10))
        pd_data = main_df["data"]
        pd_data["ks"] /= main_df["k_F"]
        if pd_data["Delta_Coulomb"][0] > 0:
            pd_data["Delta_Phonon"] *= -1
            pd_data["Delta_Coulomb"] *= -1
        plot_data = pd_data.query(f"ks > {1-OFFS[i]} & ks < {1+OFFS[i]}")
        axes[i].plot(plot_data["ks"] - 1, 1e3 * energy(plot_data["xis"] + 1e-3 * plot_data["Delta_Fock"], plot_data["Delta_Phonon"] + plot_data["Delta_Coulomb"]) - main_df["Delta_max"], label=f"$g={g}$")
        axes[i].set_xlabel(r"$k / k_\mathrm{F} - 1$")

axes[0].set_xticks([-0.0008, 0.0008])   
axes[1].set_xticks([-0.0006, 0.0006])      
axes[2].set_xticks([-0.0003, 0.0003])    
 
axes[0].set_ylabel(r"$E - \Delta_\mathrm{max} [\mathrm{meV}]$")
axes[-1].legend(loc="upper right", ncols=2)
axes[-1].set_ylim(-1, 4)
fig.subplots_adjust(wspace=0.2, hspace=0.1)
plt.show()