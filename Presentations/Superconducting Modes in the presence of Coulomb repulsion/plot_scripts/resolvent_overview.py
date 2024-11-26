import numpy as np
import matplotlib.pyplot as plt
import __path_appender as __ap
__ap.append()
from create_zoom import *
from get_data import load_panda, continuum_params
import continued_fraction_pandas as cf
from alpha_zip import *

##########################
#####   No Coulomb   #####
##########################

fig, axes = plt.subplots(nrows=2, ncols=2, sharex=True, figsize=(6.4, 6.4), sharey=True)
w_lin = np.linspace(-0.1e-3, 24e-3, 15000) + 1e-5j

for i, (ax, g) in enumerate(zip(axes[0], [0.4, 1])):
    pd_data = load_panda("continuum", "offset_20", "resolvents.json.gz", 
                        **continuum_params(N_k=20000, T=0, coulomb_scaling=0, screening=0, k_F=4.25, g=g, omega_D=10))
    resolvents = cf.ContinuedFraction(pd_data, ignore_first=5, ignore_last=120, messages=False)
    ax.set_ylim(0, 0.9)

    ax.plot(1e3 * w_lin.real, resolvents.spectral_density(w_lin, "phase_SC",     withTerminator=True), label="Phase", ls="-")
    ax.plot(1e3 * w_lin.real, resolvents.spectral_density(w_lin, "amplitude_SC", withTerminator=True), label="Higgs", ls="--")
    resolvents.mark_continuum(ax, 1e3, label=None)
    ax.set_xlim(1e3 * np.min(w_lin.real), 1e3 * np.max(w_lin.real))
    ax.text(0.08 - (i - 1) * 0.08, 0.88, f"(a.{i+1})", transform=ax.transAxes)

for i, (ax, g) in enumerate(zip(axes[1], [0.4, 1])):
    pd_data = load_panda("continuum", "theta_approx", "resolvents.json.gz", 
                        **continuum_params(N_k=8000, T=0, coulomb_scaling=0, screening=0, k_F=4.25, g=g, omega_D=10))
    resolvents = cf.ContinuedFraction(pd_data, ignore_first=5, ignore_last=120, messages=False)
    ax.set_ylim(0, 0.9)

    ax.plot(1e3 * w_lin.real, resolvents.spectral_density(w_lin, "phase_SC",     withTerminator=True), label="Phase", ls="-")
    ax.plot(1e3 * w_lin.real, resolvents.spectral_density(w_lin, "amplitude_SC", withTerminator=True), label="Higgs", ls="--")
    resolvents.mark_continuum(ax, 1e3, label=None)
    ax.set_xlim(1e3 * np.min(w_lin.real), 1e3 * np.max(w_lin.real))
    ax.text(0.08 - (i - 1) * 0.08, 0.88, f"(b.{i+1})", transform=ax.transAxes)

axes[0][0].legend(loc="upper right")
axes[1][0].set_xlabel(r"$\omega [\mathrm{meV}]$")
axes[1][1].set_xlabel(r"$\omega [\mathrm{meV}]$")
axes[0][0].set_ylabel(r"$\mathcal{A} (\omega) [\mathrm{eV}^{-1}]$")
axes[1][0].set_ylabel(r"$\mathcal{A} (\omega) [\mathrm{eV}^{-1}]$")
import os
fig.tight_layout()
fig.subplots_adjust(hspace=0, wspace=0)
fig.savefig(f"plots/{os.path.basename(__file__).split('.')[0]}.pdf")

##########################
#####  With Coulomb  #####
##########################

fig, axes = plt.subplots(nrows=2,sharex=True, figsize=(6.4, 6.4))
for i, (ax, g, label) in enumerate(alpha_label(axes, [1, 2])):
    pd_data = load_panda("continuum", "offset_10", "resolvents.json.gz", 
                        **continuum_params(N_k=20000, T=0, coulomb_scaling=1, screening=1e-4, k_F=4.25, g=g, omega_D=10))
    resolvents = cf.ContinuedFraction(pd_data, ignore_first=5, ignore_last=120, messages=False)
    ax.set_ylim(0, 0.9)

    ax.plot(1e3 * w_lin.real, resolvents.spectral_density(w_lin, "phase_SC",     withTerminator=True), label="Phase", ls="-")
    ax.plot(1e3 * w_lin.real, resolvents.spectral_density(w_lin, "amplitude_SC", withTerminator=True), label="Higgs", ls="--")
    resolvents.mark_continuum(ax, 1e3, label=None)
    ax.set_xlim(1e3 * np.min(w_lin.real), 1e3 * np.max(w_lin.real))
    ax.text(0.02, 0.88, label, transform=ax.transAxes)

axes[0].legend()
axes[1].set_xlabel(r"$\omega [\mathrm{meV}]$")
axes[0].set_ylabel(r"$\mathcal{A} (\omega) [\mathrm{eV}^{-1}]$")
axes[1].set_ylabel(r"$\mathcal{A} (\omega) [\mathrm{eV}^{-1}]$")

fig.tight_layout()
fig.subplots_adjust(hspace=0)
fig.savefig(f"plots/{os.path.basename(__file__).split('.')[0]}_coulomb.pdf")