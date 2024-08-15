import matplotlib.pyplot as plt
import numpy as np

import __path_appender as __ap
__ap.append()
from get_data import *

fig, ax = plt.subplots()

main_df = load_panda("continuum", "test5", "gap.json.gz", **continuum_params(0.0, 0.0, 4.25, 2.5, 10.))
pd_data = main_df["data"]
pd_data["ks"] /= main_df["k_F"]

energies = 1e-3 * np.sqrt((1e3 * pd_data["xis"] + 1e3 * pd_data["Delta_Fock"])**2 + (pd_data["Delta_Coulomb"] + pd_data["Delta_Phonon"])**2)# / main_df["E_F"]
ax.plot(pd_data["ks"], energies )
inner = int((main_df["discretization"] - main_df["inner_discretization"]) / 2)
ax.axvline(pd_data["ks"][inner], ls=":", color="grey")
ax.axvline(pd_data["ks"][inner + main_df["inner_discretization"]], ls=":", color="grey")

ax.set_xlabel(r"$k / k_\mathrm{F}$")
ax.set_ylabel(r"$E(k) / E_\mathrm{F}$")

fig.tight_layout()
plt.show()