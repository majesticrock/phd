import numpy as np
import matplotlib.pyplot as plt
import __path_appender as __ap
__ap.append()

from get_data import load_panda, continuum_params
import continued_fraction_pandas as cf
import plot_settings as ps

changer = [0.0, 0.1, 0.25, 0.5, 0.75, 1.0]
alphabet = "abcdefghijkl"
fig, axs = plt.subplots(len(changer), 1, sharex=True, sharey=True, figsize=(6.4, 8))

for i in range(len(changer)):
    pd_data = load_panda("continuum", "offset_10", "resolvents.json.gz", **continuum_params(8000, 0., changer[i], 1e-4, 4.25, 3.5, 10.))
    resolvents = cf.ContinuedFraction(pd_data, ignore_first=5, ignore_last=80)
    
    plotter = ps.CURVEFAMILY(6, axis=axs[i])
    plotter.set_individual_colors("nice")
    plotter.set_individual_linestyles(["-", "-.", "--", "-", "--", ":"])

    w_lin = np.linspace(0., 0.280, 15000, dtype=complex)
    w_lin += 1e-5j

    plotter.plot(1e3 * w_lin.real, resolvents.spectral_density(w_lin, "phase_SC", withTerminator=True), label="Phase")
    plotter.plot(1e3 * w_lin.real, resolvents.spectral_density(w_lin, "amplitude_SC", withTerminator=True), label="Higgs")
    resolvents.mark_continuum(axs[i], scale_factor=1e3)

    axs[len(changer) - 1].set_xlabel(r"$\omega [\mathrm{meV}]$")
    axs[i].text(0.64, 0.6, f"({alphabet[i]}) $\\alpha = {changer[i]}$", transform=axs[i].transAxes)
    axs[i].set_ylabel(r"$\mathcal{A} (\omega) [\mathrm{eV}^{-1}]$")

axs[0].set_ylim(-0.01, 1)

fig.tight_layout()
plt.show()