import matplotlib.pyplot as plt
import numpy as np

import __path_appender as __ap
__ap.append()
from get_data import *

def am_mu(ks, kf):
    return 0.12652559550141668 / ( 2 * kf ) * np.log((ks**2 + 4 * kf**2) / (ks**2))

def anderson_morel(g, omega_D, Ef, mu):
    denom = g - (mu / (1 + mu * np.log(Ef / omega_D)))
    return np.where(denom > 0, 2. * omega_D * np.exp(- 1. / denom), 0)

all_data = load_all("continuum/offset_10/N_k=20000/T=0.0", "gap.json.gz").query("k_F == 4.25 & omega_D == 10")
omega_debye = 10
k_F = 4.25
screening_factor = 0.4107320221286488672 * np.sqrt(k_F)

g_lin = np.linspace(0., 5., 200)

fig, ax = plt.subplots()

plot_data = all_data.query("coulomb_scaling == 0 & lambda_screening == 0").sort_values("g")
E_F = plot_data["E_F"].iloc[0]
ax.plot(plot_data["g"], plot_data["Delta_max"], "x", color="C0")
ax.plot(g_lin, anderson_morel(g_lin, omega_debye, E_F, 0.), color="C0")

plot_data = all_data.query("coulomb_scaling == 1 & lambda_screening == 1").sort_values("g")
E_F = plot_data["E_F"].iloc[0]
ax.plot(plot_data["g"], plot_data["Delta_max"], "x", color="C1")
ax.plot(g_lin, anderson_morel(g_lin, omega_debye, E_F, am_mu(1. * screening_factor, k_F)), color="C1")

plot_data = all_data.query("coulomb_scaling == 1 & lambda_screening == 0.0001").sort_values("g")
E_F = plot_data["E_F"].iloc[0]
ax.plot(plot_data["g"], plot_data["Delta_max"], "x", color="C2")
ax.plot(g_lin, anderson_morel(g_lin, omega_debye, E_F, am_mu(1e-4 * screening_factor, k_F)), color="C2")


ax.set_xlabel(r"$g$")
ax.set_ylabel(r"$\Delta_\mathrm{max}$ $[\mathrm{meV}]$")

plt.show()