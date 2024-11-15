import matplotlib.pyplot as plt
import numpy as np

import __path_appender as __ap
__ap.append()
from get_data import *
from scipy.optimize import curve_fit
import os
from uncertainties import ufloat
import uncertainties.umath as unp

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


dmax = 0.15
dmin = 0.01
approx = load_all("continuum/theta_approx/N_k=8000/T=0.0/coulomb_scaling=0.0", "gap.json.gz").query(
    f"k_F == 4.25 & omega_D == 10 & Delta_max < {dmax} & Delta_max > {dmin}"
    ).sort_values("g")
nc = load_all("continuum/offset_5/N_k=20000/T=0.0/coulomb_scaling=0.0", "gap.json.gz").query(
    f"k_F == 4.25 & omega_D == 10 & Delta_max < {dmax} & Delta_max > {dmin}"
    ).sort_values("g")
ls = load_all("continuum/offset_5/N_k=20000/T=0.0/coulomb_scaling=1.0", "gap.json.gz", condition="screening=1.0").query(
    f"k_F == 4.25 & omega_D == 10 & Delta_max < {dmax} & Delta_max > {dmin}"
    ).sort_values("g")
ss = load_all("continuum/offset_5/N_k=20000/T=0.0/coulomb_scaling=1.0", "gap.json.gz", condition="screening=0.0001").query(
    f"k_F == 4.25 & omega_D == 10 & Delta_max < {dmax} & Delta_max > {dmin}"
    ).sort_values("g")
dfs = [approx, nc, ls, ss]


fig, ax = plt.subplots()

for idx, (plot_data, label, screening) in enumerate(zip(dfs, ["Approx. g(k,k')$", "Full $g(k,k')$", r"$\lambda=1$", r"$\lambda=10^{-4}$"], [None, None, 1, 1e-4])):
    g_lin = np.linspace(plot_data["g"].min() - 0.015, plot_data["g"].max() + 0.02, 200)
    E_F = plot_data["E_F"].iloc[0]
    y_data = np.log(plot_data["Delta_max"] / (2. * plot_data["omega_D"].iloc[0]))
    ax.plot(plot_data["g"], y_data, ls="", marker=f"{idx+1}", ms=14, markeredgewidth=2, color=f"C{idx}", label=label)
    popt, pcov = curve_fit(fit_func, plot_data["g"], y_data, p0=(2.7, 1, (idx-1)**3 * 0.05))
    ax.plot(g_lin, fit_func(g_lin, *popt), color=f"C{idx}", ls="--")
    
    devi = np.sqrt(np.diag(pcov))
    alpha = unp.exp( ufloat(popt[0], devi[0]) )
#    print(label, ":")
#    print(f"AM prediction: 1.0000   ||   Fit: {alpha}")
#    print(f"AM prediction: 1.0000   ||   Fit: {ufloat(popt[1], devi[1])}")
#    print(f"AM prediction: {am_mu_star(E_F, screening)}   ||   Fit: {ufloat(popt[-1], devi[-1])}")
#    print(f"Bogo mu      : {bogo_mu(screening)}")

##### Results #####
#  Approx. $G$ :
#  AM prediction: 1.0000   ||   Fit: 0.99995+/-0.00034
#  AM prediction: 1.0000   ||   Fit: 0.99994+/-0.00011
#  AM prediction: 0.0      ||   Fit: (8+/-9)e-06
#  Bogo mu      : 0
#  Full $G$ :
#  AM prediction: 1.0000   ||   Fit: 0.7804+/-0.0009
#  AM prediction: 1.0000   ||   Fit: 0.9679+/-0.0004
#  AM prediction: 0.0      ||   Fit: 0.001537+/-0.000032
#  Bogo mu      : 0
#  $\lambda=1$ :
#  AM prediction: 1.0000                 ||   Fit: 0.707+/-0.012
#  AM prediction: 1.0000                 ||   Fit: 1.056+/-0.006
#  AM prediction: 0.046970974711624806   ||   Fit: 0.0526+/-0.0005
#  Bogo mu      : 0.06881081989309847
#  $\lambda=10^{-4}$ :
#  AM prediction: 1.0000               ||   Fit: 0.31+/-0.04
#  AM prediction: 1.0000               ||   Fit: 1.02+/-0.06
#  AM prediction: 0.1035079098250726   ||   Fit: 0.291+/-0.006
#  Bogo mu      : 0.3428623732280431


ax.set_xlabel(r"$g$")
ax.set_ylabel(r"$\ln (\Delta_\mathrm{max} / (2 \omega_\mathrm{D}))$")
ax.legend(ncols=2)
fig.tight_layout()
fig.savefig(f"plots/{os.path.basename(__file__).split('.')[0]}.pdf")
