import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
from matplotlib import colors
import __path_appender as __ap
__ap.append()
import continued_fraction_pandas as cf

BUILD_DIR = "python/continuum/build/"

class HeatmapPlotter:
    def __init__(self, data_frame_param, parameter_name, xlabel='Y-axis', zlabel=r'$A$ [$\mathrm{meV}^{-1}$]', title='Spectral functions'):
        data_frame = data_frame_param.sort_values(parameter_name).reset_index(drop=True)
        
        self.y_gap = np.linspace(0., 5., 2000)
        self.y_mev = np.linspace(0., 100., 2000)
        self.x = (data_frame[parameter_name]).to_numpy()
        self.resolvents = [cf.ContinuedFraction(pd_row, messages=False) for index, pd_row in data_frame.iterrows()]
        self.gaps = [2e-3 * gap for gap in data_frame["Delta_max"]]

        self.xlabel = xlabel
        self.zlabel = zlabel
        self.title = title
        
        self.plot()

    def plot(self, cmap='viridis'):
        self.fig, self.axes = plt.subplots(2, 2, sharex="col", sharey="row", figsize=(10, 8))
        self.fig.subplots_adjust(wspace=0, hspace=0.1)
        
        spectral_functions_higgs = np.array([res.spectral_density(gap * self.y_gap + 1e-4j, "amplitude_SC") for res, gap in zip(self.resolvents, self.gaps)]).transpose()
        spectral_functions_phase = np.array([res.spectral_density(gap * self.y_gap + 1e-4j, "phase_SC") for res, gap in zip(self.resolvents, self.gaps)]).transpose()
        spectral_functions_higgs_mev = np.array([res.spectral_density(1e-3 * self.y_mev + 1e-4j, "amplitude_SC") for res in self.resolvents]).transpose()
        spectral_functions_phase_mev = np.array([res.spectral_density(1e-3 * self.y_mev + 1e-4j, "phase_SC") for res in self.resolvents]).transpose()
        
        norm = colors.Normalize(vmin=0, vmax=1)
        
        contour = self.axes[0][0].pcolormesh(self.x, self.y_gap, spectral_functions_higgs, cmap=cmap, norm=norm)
        self.axes[0][1].pcolormesh(self.x, self.y_gap, spectral_functions_phase, cmap=cmap, norm=norm)
        self.axes[1][0].pcolormesh(self.x, self.y_mev, spectral_functions_higgs_mev, cmap=cmap, norm=norm)
        self.axes[1][1].pcolormesh(self.x, self.y_mev, spectral_functions_phase_mev, cmap=cmap, norm=norm)
        
        cbar = self.fig.colorbar(contour, ax=self.axes, orientation='vertical', fraction=0.046, pad=0.04, extend='max')
        cbar.set_label(self.zlabel)

        self.axes[0][0].set_ylabel(r"$\omega$ $[2 \Delta]$")
        self.axes[1][0].set_ylabel(r"$\omega$ [meV]")
        
        self.axes[1][0].set_xlabel(self.xlabel)
        self.axes[1][1].set_xlabel(self.xlabel)
        
        self.axes[0][0].set_title(r"$\mathcal{A}_\mathrm{Higgs}$")
        self.axes[0][1].set_title(r"$\mathcal{A}_\mathrm{Phase}$")
        self.fig.suptitle(self.title)

    def show_plot(self):
        plt.show()

    def save_plot(self, filename):
        self.fig.savefig(filename)


import __path_appender as __ap
__ap.append()
from get_data import *

all_data = load_all("continuum/offset_10/N_k=20000/T=0.0", "resolvents.json.gz")


HeatmapPlotter(all_data.query("coulomb_scaling == 1 & k_F == 4.25 & omega_D == 10 & g == 0.5"), "lambda_screening", xlabel=r"$\lambda$", title=r"$g = 0.5$")



omega_D_small_screening = HeatmapPlotter(all_data.query("coulomb_scaling == 1 & lambda_screening == 0.0001 & k_F == 4.25 & g == 1 & omega_D < 21"), "omega_D", xlabel=r"$\omega_\mathrm{D}$ [meV]", title=r"$\lambda = 0.0001$")
omega_D_large_screening = HeatmapPlotter(all_data.query("coulomb_scaling == 1 & lambda_screening == 1      & k_F == 4.25 & g == 1 & omega_D < 21"), "omega_D", xlabel=r"$\omega_\mathrm{D}$ [meV]", title=r"$\lambda = 1$")
omega_D_no_coulomb      = HeatmapPlotter(all_data.query("coulomb_scaling == 0 & lambda_screening == 0      & k_F == 4.25 & g == 1 & omega_D < 21"), "omega_D", xlabel=r"$\omega_\mathrm{D}$ [meV]", title="No Coulomb")
omega_D_small_screening.save_plot(BUILD_DIR + "omega_D_small_screening.pdf")
omega_D_large_screening.save_plot(BUILD_DIR + "omega_D_large_screening.pdf")
omega_D_no_coulomb.save_plot(     BUILD_DIR + "omega_D_no_coulomb.pdf")

g_small_screening = HeatmapPlotter(all_data.query("coulomb_scaling == 1 & lambda_screening == 0.0001 & k_F == 4.25 & omega_D == 10 & g > 0.7 & g < 5"),  "g", xlabel=r"$g$", title=r"$\lambda = 0.0001$")
g_large_screening = HeatmapPlotter(all_data.query("coulomb_scaling == 1 & lambda_screening == 1      & k_F == 4.25 & omega_D == 10 & g > 0.7 & g < 5"),  "g", xlabel=r"$g$", title=r"$\lambda = 1$")
g_no_coulomb      = HeatmapPlotter(all_data.query("coulomb_scaling == 0 & lambda_screening == 0      & k_F == 4.25 & omega_D == 10 & g > 0.25 & g < 4"), "g", xlabel=r"$g$", title="No Coulomb")
g_small_screening.save_plot(BUILD_DIR + "g_small_screening.pdf")
g_large_screening.save_plot(BUILD_DIR + "g_large_screening.pdf")
g_no_coulomb.save_plot(     BUILD_DIR + "g_no_coulomb.pdf")

plt.show()