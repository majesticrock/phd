import numpy as np
import matplotlib.pyplot as plt
import __path_appender as __ap
__ap.append()
from create_zoom import *
from get_data import load_panda, continuum_params
pd_data = load_panda("continuum", "test", "resolvents.json.gz",
                    **continuum_params(N_k=4000, T=0, coulomb_scaling=0, screening=1e-4, k_F=4.25, g=0.5, omega_D=10))

import continued_fraction_pandas as cf
import plot_settings as ps

resolvents = cf.ContinuedFraction(pd_data, ignore_first=30, ignore_last=90)

fig, ax = plt.subplots()
ax.set_ylim(-0.05, 1)
ax.set_xlabel(r"$\omega [\mathrm{meV}]$")
ax.set_ylabel(r"$\mathcal{A} (\omega) [\mathrm{eV}^{-1}]$")

plotter = ps.CURVEFAMILY(6, axis=ax)
plotter.set_individual_colors("nice")
plotter.set_individual_linestyles(["-", "-.", "--", "-", "--", ":"])

w_lin = np.linspace(-0.005 * pd_data["continuum_boundaries"][1], 1.1 * pd_data["continuum_boundaries"][1], 150000, dtype=complex)
#w_lin = np.linspace(0, 150, 15000, dtype=complex)
w_lin += 1e-4j

plotter.plot(1e3 * w_lin.real, resolvents.spectral_density(w_lin, "phase_SC",     withTerminator=True), label="Phase")
plotter.plot(1e3 * w_lin.real, resolvents.spectral_density(w_lin, "amplitude_SC", withTerminator=True), label="Higgs")

resolvents.mark_continuum(ax, 1e3)

#import gzip
#with gzip.open("data/continuum/test/full_diag/-values.dat.gz", 'rt') as f_open:
#    M_ev = np.loadtxt(f_open)
#with gzip.open("data/continuum/test/full_diag/sc-weights.dat.gz", 'rt') as f_open:
#    M_w = np.loadtxt(f_open)
#    
#def resolvent(z, evs, weights):
#    ret = np.zeros(len(z), dtype=complex)
#    for i in range(len(evs)):
#        ret += weights[i] / (z - evs[i])
#    return - ret.imag / np.pi
#plotter.plot(1e3 * w_lin.real, resolvent(w_lin, M_ev, M_w), label="Exact")

#axins = create_zoom(ax, 0.2, 0.3, 0.4, 0.4, xlim=(0., 3e4 * pd_data["continuum_boundaries"][0]), ylim=(0., 0.08))
#resolvents.mark_continuum(axins, 1e3)

ax.set_xlim(1e3 * np.min(w_lin.real), 1e3 * np.max(w_lin.real))
ax.legend()
fig.tight_layout()
plt.show()