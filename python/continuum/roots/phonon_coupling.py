import enum
import matplotlib.pyplot as plt
import numpy as np

import __path_appender as __ap
__ap.append()
from get_data import *
from scipy.interpolate import InterpolatedUnivariateSpline
from ez_fit import *

TYPE = "g"
LABEL= r"$g$"

pd_data = load_all("continuum/exact_2000", "gap.json.gz").query('omega_D == 10 & E_F == 9.3 & \
                    (coulomb_scaling == 1 | coulomb_scaling == 0.75 | coulomb_scaling == 0.5 | coulomb_scaling == 0.25)')
pd_data.sort_values(TYPE, inplace=True)
pd_data.reset_index(inplace=True)

coulomb_group = pd_data.groupby("coulomb_scaling")

min_positions = np.empty(coulomb_group.ngroups, dtype=np.ndarray)
min_energies  = np.empty(coulomb_group.ngroups, dtype=np.ndarray)

fig, ax = plt.subplots()

for coulomb_idx, (coulomb_scaling, iter_frame) in enumerate(coulomb_group):
    DATA_SIZE = len(iter_frame.reset_index().index)
    min_positions[coulomb_idx] = np.zeros(DATA_SIZE)
    min_energies[coulomb_idx]  = np.zeros(DATA_SIZE)
    
    for index, pd_row in iter_frame.reset_index().iterrows():
        data_arrays = pd_row["data"]
        HALF_LEN = int(len(data_arrays["ks"]) / 2)
        data_arrays["ks"] -= pd_row["k_F"]
        delta_arr = data_arrays["Delta_Phonon"][:HALF_LEN] + data_arrays["Delta_Coulomb"][:HALF_LEN]

        f = InterpolatedUnivariateSpline( data_arrays["ks"][:HALF_LEN], delta_arr, k=3)
        cr_pts = f.roots()
        min_index = np.argmin(cr_pts)
        min_positions[coulomb_idx][index] = cr_pts[min_index]

        #f = InterpolatedUnivariateSpline( data_arrays["xis"][:HALF_LEN], delta_arr, k=3)
        #cr_pts = f.roots()
        #min_index = np.argmin(cr_pts)
        #min_energies[coulomb_idx][index] = cr_pts[min_index]

    ax.plot(iter_frame[TYPE],       min_positions[coulomb_idx] / np.sqrt(1e-3 * pd_row["Delta_max"]), "v", c=f"C{coulomb_idx}", label=fr"$\alpha = {coulomb_scaling}$")
    ez_linear_fit(iter_frame[TYPE], min_positions[coulomb_idx] / np.sqrt(1e-3 * pd_row["Delta_max"]), ax, c=f"C{coulomb_idx}")

ax.set_ylabel(r"$k_0 - k_\mathrm{F} [\sqrt{\mathrm{eV}}]$")
ax.legend()
ax.set_xlabel(f"{LABEL} [eV]")
fig.tight_layout()

import os
plt.savefig(f"python/build/root_{os.path.basename(__file__).split('.')[0]}.svg")
plt.show()