import matplotlib.pyplot as plt
import numpy as np

import __path_appender as __ap
__ap.append()
from get_data import *
from scipy.optimize import curve_fit
import os

k_F = 4.25
screening_factor = 0.4107320221286488672 * np.sqrt(k_F)
omega_D = 0.01

def bogo_mu(lam):
    if lam is None or lam == 0:
        return 0
    ks = lam * screening_factor
    #0.12652559550141668 sqrt(eV) = 1 / (4 * pi * pi * epsilon_0)
    return 0.12652559550141668 / ( 2 * k_F ) * np.log((ks**2 + 4 * k_F**2) / (ks**2))

def am_mu_star(Ef, lam):
    mu = bogo_mu(lam)
    return mu / (1 + mu * np.log(Ef / omega_D))

def anderson_morel(g, Ef, mu):
    denom = g - am_mu_star(Ef, mu)
    return np.where(denom > 0, 2. * omega_D * np.exp(- 1. / denom), 0)

def fit_func(g, lnA, B, mu_star):
    return lnA - B / (g - mu_star)


dmax = 0.3
nc = load_all("continuum/offset_5/N_k=20000/T=0.0/coulomb_scaling=0.0", "gap.json.gz").query(
    f"k_F == 4.25 & omega_D == 10 & Delta_max < {dmax} & Delta_max > 0.008"
    ).sort_values("g")
ls = load_all("continuum/offset_5/N_k=20000/T=0.0/coulomb_scaling=1.0", "gap.json.gz", condition="screening=1.0").query(
    f"k_F == 4.25 & omega_D == 10 & Delta_max < {dmax} & Delta_max > 0.008"
    ).sort_values("g")
ss = load_all("continuum/offset_5/N_k=20000/T=0.0/coulomb_scaling=1.0", "gap.json.gz", condition="screening=0.0001").query(
    f"k_F == 4.25 & omega_D == 10 & Delta_max < {dmax} & Delta_max > 0.008"
    ).sort_values("g")
dfs = [nc, ls, ss]


fig, ax = plt.subplots()

for idx, (plot_data, label, screening) in enumerate(zip(dfs, ["No Coulomb", r"$\lambda=1$", r"$\lambda=10^{-4}$"], [None, 1, 1e-4])):
    g_lin = np.linspace(plot_data["g"].min() - 0.02, plot_data["g"].max() + 0.02, 200)
    E_F = plot_data["E_F"].iloc[0]
    y_data = np.log(plot_data["Delta_max"] / (2. * plot_data["omega_D"].iloc[0]))
    ax.plot(plot_data["g"], y_data, ls="", marker=f"{idx+1}", ms=12, color=f"C{idx}", label=label)
    popt, pcov = curve_fit(fit_func, plot_data["g"], y_data, p0=(2.7, 1, idx**3 * 0.05))
    ax.plot(g_lin, fit_func(g_lin, *popt), color=f"C{idx}")
    
    #devi = np.sqrt(np.diag(pcov))
    #print(label, ":")
    #print(f"AM prediction: 0.0000   ||   Fit: {popt[0]:1.4f} +/- {devi[0]:1.4f}")
    #print(f"AM prediction: 1.0000   ||   Fit: {popt[1]:1.4f} +/- {devi[1]:1.4f}")
    #print(f"AM prediction: {am_mu_star(E_F, screening):1.4f}   ||   Fit: {popt[2]:1.4f} +/- {devi[2]:1.4f}")

##### Results #####
#  No Coulomb :
#  AM prediction: 0.0000   ||   Fit: -0.2256 +/- 0.0099
#  AM prediction: 1.0000   ||   Fit: 0.9766 +/- 0.0035
#  AM prediction: 0.0000   ||   Fit: 0.0007 +/- 0.0003
#  $\lambda=1$ :
#  AM prediction: 0.0000   ||   Fit: -0.3041 +/- 0.0072
#  AM prediction: 1.0000   ||   Fit: 1.0715 +/- 0.0029
#  AM prediction: 0.0470   ||   Fit: 0.0514 +/- 0.0003
#  $\lambda=10^{-4}$ :
#  AM prediction: 0.0000   ||   Fit: -1.0022 +/- 0.1147
#  AM prediction: 1.0000   ||   Fit: 1.0651 +/- 0.0524
#  AM prediction: 0.1035   ||   Fit: 0.2885 +/- 0.0054


ax.set_xlabel(r"$g$")
ax.set_ylabel(r"$\ln (\Delta_\mathrm{max} / (2 \omega_\mathrm{D}))$")
ax.legend()
fig.tight_layout()
fig.savefig(f"plots/{os.path.basename(__file__).split('.')[0]}.pdf")
