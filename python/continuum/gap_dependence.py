import matplotlib.pyplot as plt
import numpy as np

import __path_appender as __ap
__ap.append()
from get_data import *

k_F = 4.25
screening_factor = 0.4107320221286488672 * np.sqrt(k_F)
omega_D = 0.01

def am_mu(ks):
    return 0.12652559550141668 / ( 2 * k_F ) * np.log((ks**2 + 4 * k_F**2) / (ks**2))

def mu_star(Ef, mu):
    return mu / (1 + mu * np.log(Ef / omega_D))

def anderson_morel(g, Ef, mu):
    denom = g - mu_star(Ef, mu)
    return np.where(denom > 0, 2. * omega_D * np.exp(- 1. / denom), 0)

def x_func(g, offset):
    return (g - offset)

ss = load_all("continuum/offset_5/N_k=8000/T=0.0/coulomb_scaling=1.0", "gap.json.gz", condition="screening=0.0001").query(
    "k_F == 4.25 & omega_D == 10 & Delta_max < 0.4 & Delta_max > 0"
    )
ls = load_all("continuum/offset_5/N_k=8000/T=0.0/coulomb_scaling=1.0", "gap.json.gz", condition="screening=1.0").query(
    "k_F == 4.25 & omega_D == 10 & Delta_max < 0.4 & Delta_max > 0"
    )
nc = load_all("continuum/offset_5/N_k=8000/T=0.0/coulomb_scaling=0.0", "gap.json.gz").query(
    "k_F == 4.25 & omega_D == 10 & Delta_max < 0.4 & Delta_max > 0"
    )

g_lin = np.linspace(0., 0.8, 200)

fig, ax = plt.subplots()

plot_data = nc.query("coulomb_scaling == 0 & lambda_screening == 0").sort_values("g")
E_F = plot_data["E_F"].iloc[0]
ax.plot(x_func(plot_data["g"], 0), np.log(plot_data["Delta_max"]), "x", color="C0")
#ax.plot(g_lin, anderson_morel(g_lin, E_F, 0.), color="C0")

plot_data = ls.query("coulomb_scaling == 1 & lambda_screening == 1").sort_values("g")
E_F = plot_data["E_F"].iloc[0]
ax.plot(x_func(plot_data["g"], 0.09), np.log(plot_data["Delta_max"]), "x", color="C1")
#ax.plot(g_lin, anderson_morel(g_lin, E_F, am_mu(1. * screening_factor)), color="C1")

plot_data = ss.query("coulomb_scaling == 1 & lambda_screening == 0.0001").sort_values("g")
E_F = plot_data["E_F"].iloc[0]
ax.plot(x_func(plot_data["g"], 0.375), np.log(plot_data["Delta_max"]), "x", color="C2")
#ax.plot(g_lin, anderson_morel(g_lin, E_F, am_mu(1e-4 * screening_factor)), color="C2")

ax.set_xlabel(r"$g$")
ax.set_ylabel(r"$\Delta_\mathrm{max}$ $[\mathrm{meV}]$")

plt.show()