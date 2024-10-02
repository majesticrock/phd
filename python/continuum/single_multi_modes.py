import numpy as np
import matplotlib.pyplot as plt
import __path_appender as __ap
__ap.append()

from get_data import load_panda, continuum_params
import continued_fraction_pandas as cf
import plot_settings as ps

screenings = [1e-4, 1e-2, 0.1, 0.25, 0.5, 0.75, 1., 1.5, 10., 50., 100., r"\lambda"]
phonons = [0.6, 0.8, 1., 1.2, 1.4, 2., "g"]
kfs = [3., 3.5, 4.0, 4.5, 5.0, "k_\mathrm{F}"]
Ts = [11., 12., 13., 14., "T"]

plotted = screenings
n_plots = len(plotted) - 1

alphabet = "abcdefghijkl"
fig, axs = plt.subplots(n_plots, 1, sharex=True, sharey=True, figsize=(6.4, 12))

for i in range(n_plots):
    pd_data = load_panda("continuum", "offset_10", "resolvents.json.gz", 
                #**continuum_params(N_k=8000, T=0.0, coulomb_scaling=1., screening=1e-4, k_F=4.25, g=phonons[i], omega_D=10.)
                **continuum_params(N_k=8000, T=0.0, coulomb_scaling=1., screening=screenings[i], k_F=4.25, g=0.6, omega_D=10.)
                #**continuum_params(N_k=8000, T=0.0, coulomb_scaling=1., screening=1e-4, k_F=kfs[i], g=1., omega_D=10.)
                #**continuum_params(N_k=8000, T=0.0, coulomb_scaling=0., screening=1e-4, k_F=kfs[i], g=0.5, omega_D=10.)
                #**continuum_params(N_k=8000, T=Ts[i], coulomb_scaling=1., screening=1e-4, k_F=4.25, g=1., omega_D=10.)
                )
    resolvents = cf.ContinuedFraction(pd_data, ignore_first=5, ignore_last=80)
    
    plotter = ps.CURVEFAMILY(6, axis=axs[i])
    plotter.set_individual_colors("nice")
    plotter.set_individual_linestyles(["-", "-.", "--", "-", "--", ":"])

    w_lin = np.linspace(0., 0.080, 15000, dtype=complex)
    w_lin += 1e-5j

    plotter.plot(1e3 * w_lin.real, resolvents.spectral_density(w_lin, "phase_SC", withTerminator=True), label="Phase")
    plotter.plot(1e3 * w_lin.real, resolvents.spectral_density(w_lin, "amplitude_SC", withTerminator=True), label="Higgs")
    resolvents.mark_continuum(axs[i], scale_factor=1e3)

    axs[n_plots - 1].set_xlabel(r"$\omega [\mathrm{meV}]$")
    axs[i].text(0.64, 0.6, f"({alphabet[i]}) ${plotted[-1]} = {plotted[i]}$", transform=axs[i].transAxes)
    axs[i].set_ylabel(r"$\mathcal{A} (\omega) [\mathrm{eV}^{-1}]$")

axs[0].set_ylim(-0.01, 1)
axs[0].set_xlim(1e3 * w_lin.real.min(), 1e3 * w_lin.real.max())

fig.tight_layout()
plt.show()