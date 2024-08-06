import numpy as np
import matplotlib.pyplot as plt
import __path_appender as __ap
__ap.append()

from get_data import load_panda, continuum_params
pd_data = load_panda("continuum/", "test", "resolvents.json.gz", **continuum_params(0., 1., 1.3, 5., 10.))

import continued_fraction_pandas as cf
import plot_settings as ps

resolvents = cf.ContinuedFraction(pd_data, ignore_first=5, ignore_last=60)

fig, ax = plt.subplots()
#ax.set_ylim(-0.05, 1e3)
ax.set_xlabel(r"$\omega [\mathrm{eV}]$")
ax.set_ylabel(r"$\mathcal{A} (\omega) [\mathrm{eV}^{-1}]$")

plotter = ps.CURVEFAMILY(6, axis=ax)
plotter.set_individual_colors("nice")
plotter.set_individual_linestyles(["-", "-.", "--", "-", "--", ":"])

w_lin = np.linspace(-0.01 * pd_data["continuum_boundaries"][0], 1.1 * pd_data["continuum_boundaries"][1], 5000, dtype=complex)
w_lin += 1e-6j

plotter.plot(w_lin.real, resolvents.spectral_density(w_lin, "phase_SC"), label="Phase")
plotter.plot(w_lin.real, resolvents.spectral_density(w_lin, "amplitude_SC"), label="Higgs")

resolvents.mark_continuum(ax)

ax.legend()
fig.tight_layout()
plt.show()