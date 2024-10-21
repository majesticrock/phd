import numpy as np
import matplotlib.pyplot as plt
import __path_appender as __ap
__ap.append()

from get_data import load_panda, hubbard_params
pd_data = load_panda("hubbard/square", "test", "resolvents.json.gz", **hubbard_params(0.0, -2.5, 0.0))

import continued_fraction_pandas as cf
import plot_settings as ps

resolvents = cf.ContinuedFraction(pd_data)

fig, ax = plt.subplots()
ax.set_ylim(-0.05, 1.)
ax.set_xlabel(r"$\omega [t]$")
ax.set_ylabel(r"$\mathcal{A} (\omega) [t^{-1}]$")

plotter = ps.CURVEFAMILY(6, axis=ax)
plotter.set_individual_colors("nice")
plotter.set_individual_linestyles(["-", "-.", "--", "-", "--", ":"])

w_lin = np.linspace(-0.01, pd_data["continuum_boundaries"][1] + 0.3, 5000, dtype=complex)
w_lin += 1e-6j

plotter.plot(w_lin, resolvents.spectral_density(w_lin, "phase_SC_a"), label="Phase")
plotter.plot(w_lin, resolvents.spectral_density(w_lin, "amplitude_SC_a"), label="Higgs")
plotter.plot(w_lin, resolvents.spectral_density(w_lin, "amplitude_CDW_a"), label="CDW")
plotter.plot(w_lin, resolvents.spectral_density(w_lin, "amplitude_AFM_a"), label="l.AFM")
plotter.plot(w_lin, resolvents.spectral_density(w_lin, "amplitude_AFM_transversal_a"), label="t.AFM")

resolvents.mark_continuum(ax)

ax.legend()
fig.tight_layout()
plt.show()