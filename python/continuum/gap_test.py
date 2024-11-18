import matplotlib.pyplot as plt
import numpy as np

import __path_appender as __ap
__ap.append()
from get_data import *
from create_zoom import *

X_BOUNDS = [-0.1, 0.1]

fig, ax = plt.subplots()

main_df = load_panda("continuum", "offset_10", "gap.json.gz",
                    **continuum_params(N_k=20000, T=0, coulomb_scaling=0, screening=1e-4, k_F=4.25, g=0.6, omega_D=10))
pd_data = main_df["data"]
pd_data["ks"] /= main_df["k_F"]
pd_data["xis"] /= main_df["E_F"]

if "imag_Delta_Phonon" in pd_data:
    phonon = pd_data["Delta_Phonon"].to_numpy() + 1j * pd_data["imag_Delta_Phonon"].to_numpy()
    coulomb = pd_data["Delta_Coulomb"].to_numpy() + 1j * pd_data["imag_Delta_Coulomb"].to_numpy()
    ax.plot(pd_data["ks"], np.real(phonon + coulomb), c="blue")
    ax2 = ax.twinx()
    ax2.plot(pd_data["ks"], np.imag(phonon + coulomb), c="red", ls="--")
    
    ax.tick_params(axis='y', colors='blue')
    ax2.tick_params(axis='y', colors='red')
else:
    if pd_data["Delta_Coulomb"][0] > 0:
        pd_data["Delta_Phonon"] *= -1
        pd_data["Delta_Coulomb"] *= -1
    ax.plot(pd_data["ks"], pd_data["Delta_Phonon"] + pd_data["Delta_Coulomb"], "k-", label=r"$\Delta$")
    pd_data.plot(x="ks", y=["Delta_Phonon", "Delta_Coulomb", "Delta_Fock"], ax=ax, style=['--', '--', ':'], label=[r"$\Delta_\mathrm{Phonon}$", r"$\Delta_\mathrm{Coulomb}$", r"$\epsilon_\mathrm{C}$"])

inner = int((main_df["discretization"] - main_df["inner_discretization"]) / 2)

#axins = create_zoom(ax, 0.1, 0.35, 0.3, 0.59, xlim=(1-0.005, 1.005), ylim=(1.2 * np.min(pd_data["Delta_Coulomb"]), 1.05 * np.max(pd_data["Delta_Phonon"])))

ax.set_xlabel(r"$k / k_\mathrm{F}$")
ax.set_ylabel(r"$\Delta [\mathrm{meV}]$")

ax.legend()
fig.tight_layout()

plt.show()