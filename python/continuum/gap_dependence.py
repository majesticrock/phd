import matplotlib.pyplot as plt
import numpy as np

import __path_appender as __ap
__ap.append()
from get_data import *
from ez_fit import * 

k_F = 4.25
screening_factor = 0.4107320221286488672 * np.sqrt(k_F)
omega_D = 0.01

def am_mu(ks):
    #0.12652559550141668 sqrt(eV) = 1 / (4 * pi * pi * epsilon_0)
    return 0.12652559550141668 / ( 2 * k_F ) * np.log((ks**2 + 4 * k_F**2) / (ks**2))

def mu_star(Ef, mu):
    return mu / (1 + mu * np.log(Ef / omega_D))

def anderson_morel(g, Ef, mu):
    denom = g - mu_star(Ef, mu)
    return np.where(denom > 0, 2. * omega_D * np.exp(- 1. / denom), 0)

def x_func(g, offset):
    return  1/ (g - offset)

dmax = 0.4

ss = load_all("continuum/offset_5/N_k=20000/T=0.0/coulomb_scaling=1.0", "gap.json.gz", condition="screening=0.0001").query(
    f"k_F == 4.25 & omega_D == 10 & Delta_max < {dmax} & Delta_max > 0.008"
    )
ls = load_all("continuum/offset_5/N_k=20000/T=0.0/coulomb_scaling=1.0", "gap.json.gz", condition="screening=1.0").query(
    f"k_F == 4.25 & omega_D == 10 & Delta_max < {dmax} & Delta_max > 0.008"
    )
nc = load_all("continuum/offset_5/N_k=20000/T=0.0/coulomb_scaling=0.0", "gap.json.gz").query(
    f"k_F == 4.25 & omega_D == 10 & Delta_max < {dmax} & Delta_max > 0.008"
    )

g_lin = np.linspace(0., 0.8, 200)

fig, ax = plt.subplots()

plot_data = nc.query("coulomb_scaling == 0 & lambda_screening == 0").sort_values("g")
E_F = plot_data["E_F"].iloc[0]
x_data = x_func(plot_data["g"], 0)
ax.plot(x_data, np.log(plot_data["Delta_max"]), "X", color="C0", label="No Coulomb")
popt, pcov = ez_linear_fit(x_data, np.log(plot_data["Delta_max"]), ax, 
              x_bounds=[0.8 * min(x_data), 1.2 * max(x_data)], color="C0")
print(popt)




plot_data = ls.query("coulomb_scaling == 1 & lambda_screening == 1").sort_values("g")
E_F = plot_data["E_F"].iloc[0]
x_data = x_func(plot_data["g"], mu_star(E_F, am_mu(screening_factor)))
ax.plot(x_data, np.log(plot_data["Delta_max"]), "X", color="C1", label=r"$\lambda=1$")
popt, pcov = ez_linear_fit(x_data, np.log(plot_data["Delta_max"]), ax, 
              x_bounds=[0.8 * min(x_data), 1.2 * max(x_data)], color="C1")
print(popt)



plot_data = ss.query("coulomb_scaling == 1 & lambda_screening == 0.0001").sort_values("g")
E_F = plot_data["E_F"].iloc[0]
x_data = x_func(plot_data["g"], 2.85 * mu_star(E_F, am_mu(1e-4 * screening_factor)))
ax.plot(x_data, np.log(plot_data["Delta_max"]), "X", color="C2", label=r"$\lambda=10^{-4}$")
popt, pcov = ez_linear_fit(x_data, np.log(plot_data["Delta_max"]), ax, 
              x_bounds=[0.8 * min(x_data), 1.2 * max(x_data)], color="C2")
print(popt)



ax.set_xlabel(r"$1 / (g - \mu^*)$")
ax.set_ylabel(r"$\ln (\Delta_\mathrm{max}$ / $\mathrm{meV})$")
ax.legend()
fig.tight_layout()

plt.show()