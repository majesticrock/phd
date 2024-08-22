import numpy as np
import matplotlib.pyplot as plt
import __path_appender as __ap
__ap.append()

from get_data import load_panda, continuum_params
import continued_fraction_pandas as cf
import plot_settings as ps

N_ks = [2000, 4000, 8000]
couplings = [1.5, 2.5, 3.5, 4.5]
fig, axs = plt.subplots(len(N_ks), len(couplings), sharex="col", sharey=True, figsize=(9, 6))

for i in range(len(N_ks)):
    for j in range(len(couplings)):
        pd_data = load_panda("continuum", "offset_20", "resolvents.json.gz", **continuum_params(N_ks[i], 0., 0., 1e-4, 4.25, couplings[j], 10.))
        resolvents = cf.ContinuedFraction(pd_data, ignore_first=5, ignore_last=80)
        
        plotter = ps.CURVEFAMILY(6, axis=axs[i][j])
        plotter.set_individual_colors("nice")
        plotter.set_individual_linestyles(["-", "-.", "--", "-", "--", ":"])

        w_lin = np.linspace(-0.01 * pd_data["continuum_boundaries"][0], 1.1 * pd_data["continuum_boundaries"][1], 15000, dtype=complex)
        w_lin += 1e-5j

        plotter.plot(1e3 * w_lin.real, resolvents.spectral_density(w_lin, "phase_SC", withTerminator=True), label="Phase")
        plotter.plot(1e3 * w_lin.real, resolvents.spectral_density(w_lin, "amplitude_SC", withTerminator=True), label="Higgs")

        resolvents.mark_continuum(axs[i][j], scale_factor=1e3)

        axs[len(N_ks) - 1][j].set_xlabel(r"$\omega [\mathrm{meV}]$")
        axs[0][j].set_xlim(0., 1.1e3 * pd_data["continuum_boundaries"][0])
        axs[i][j].text(0.2, 0.85, f"$g={couplings[j]} \\mathrm{{meV}}$", transform=axs[i][j].transAxes)
        axs[i][j].text(0.2, 0.7, f"$N_k={int(N_ks[i] / 2)}$", transform=axs[i][j].transAxes)
    axs[i][0].set_ylabel(r"$\mathcal{A} (\omega) [\mathrm{eV}^{-1}]$")
        
axs[0][0].set_ylim(-0.01, 5)

fig.tight_layout()
plt.show()