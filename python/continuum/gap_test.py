import matplotlib.pyplot as plt
import numpy as np

import __path_appender as __ap
__ap.append()
from create_zoom import create_zoom
from get_data import *

X_BOUNDS = [-0.1, 0.1]

fig, ax = plt.subplots()

main_df = load_panda("continuum", "test", "gap.json.gz", **continuum_params(0.0, 1.0, 0.45, 5., 10.))
pd_data = main_df["data"]

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
    ax.plot(pd_data["ks"], pd_data["Delta_Phonon"] + pd_data["Delta_Coulomb"], "k-", label=r"$\Delta_\mathrm{SC}$")
    pd_data.plot(x="ks", y=["Delta_Phonon", "Delta_Coulomb", "Delta_Fock"], ax=ax, style=['--', '--', '-'], label=[r"$\Delta_\mathrm{Phonon}$", r"$\Delta_\mathrm{Coulomb}$", r"$\Delta_\mathrm{Fock}$"])

inner = int((main_df["discretization"] - main_df["inner_discretization"]) / 2)
ax.axvline(pd_data["ks"][inner], ls=":", color="grey")
ax.axvline(pd_data["ks"][inner + main_df["inner_discretization"]], ls=":", color="grey")
#(pd_data["Delta_Phonon"] + pd_data["Delta_Coulomb"]).plot(ax=ax, x="ks", style="k-", label=r"$\Delta_\mathrm{SC}$")

#from scipy import integrate
#ax.text(0.4, 0.15, fr'$\int \Delta_\mathrm{{Phonon}} \mathrm{{d}}k / k_\mathrm{{F}} = {integrate.trapezoid(pd_data["Delta_Phonon"], x=pd_data["ks"]) / main_df["k_F"]:.5f}$ meV', transform=ax.transAxes)
#ax.text(0.4, 0.05, fr'$\int \Delta_\mathrm{{Coulomb}} \mathrm{{d}}k / k_\mathrm{{F}} = {integrate.trapezoid(pd_data["Delta_Coulomb"], x=pd_data["ks"]) / main_df["k_F"]:.5f}$ meV', transform=ax.transAxes)
#axins = create_zoom(ax, 0.4, 0.3, 0.3, 0.6, xlim=(1-0.02, 1.02), ylim=(1.1 * np.min(pd_data["Delta_Coulomb"]), 1.05 * np.max(pd_data["Delta_Phonon"])))

ax.set_xlabel(r"$k / k_\mathrm{F}$")
ax.set_ylabel(r"$\Delta [\mathrm{meV}]$")

ax.legend()
fig.tight_layout()

import os
plt.savefig(f"python/build/{os.path.basename(__file__).split('.')[0]}.svg")
plt.show()