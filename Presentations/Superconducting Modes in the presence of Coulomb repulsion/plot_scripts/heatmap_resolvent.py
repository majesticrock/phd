import matplotlib.pyplot as plt
import numpy as np
import __path_appender as __ap
__ap.append()
import continued_fraction_pandas as cf
import os
from get_data import *

BUILD_DIR = "plots/"
FILE_ENDING = ".pdf"

class HeatmapPlotter:
    def __init__(self, data_frame_param, parameter_name, xlabel='Y-axis', zlabel=r'$A$ [$\mathrm{meV}^{-1}$]', xscale="linear", yscale="linear"):
        data_frame = data_frame_param.sort_values(parameter_name).reset_index(drop=True)
        
        self.y = np.linspace(0., 60., 2000)
        self.x = (data_frame[parameter_name]).to_numpy()
        self.resolvents = [cf.ContinuedFraction(pd_row, messages=False) for index, pd_row in data_frame.iterrows()]
        self.gaps = [2 * gap for gap in data_frame["Delta_max"]]

        self.xlabel = xlabel
        self.zlabel = zlabel
        self.xscale = xscale
        self.yscale = yscale

    def plot(self, axes, cmap='inferno', labels=True):
        spectral_functions_higgs = np.array([res.spectral_density(1e-3 * self.y + 1e-6j, "amplitude_SC") for res in self.resolvents]).transpose()
        spectral_functions_phase = np.array([res.spectral_density(1e-3 * self.y + 1e-6j, "phase_SC") for res in self.resolvents]).transpose()

        vmax = max(spectral_functions_higgs.max(), spectral_functions_phase.max())
        levels = np.linspace(0., min(1., vmax), 101, endpoint=True)

        contour_higgs = axes[0].contourf(self.x, self.y, spectral_functions_higgs, cmap=cmap, levels=levels, extend='max', zorder=-20)
        contour_phase = axes[1].contourf(self.x, self.y, spectral_functions_phase, cmap=cmap, levels=levels, extend='max', zorder=-20)
        contour_higgs.set_edgecolor('face')
        contour_phase.set_edgecolor('face')
        
        for ax in axes:
            ax.plot(self.x, self.gaps, color="lime", ls=":")
            ax.set_rasterization_zorder(-10)
            ax.set_ylim(0., max(self.y))
            ax.set_xscale(self.xscale)
            ax.set_yscale(self.yscale)

        if(labels):
            axes[0].set_ylabel(r"$\omega [\mathrm{eV}]$")
            axes[1].set_ylabel(r"$\omega [\mathrm{eV}]$")
        axes[1].set_xlabel(self.xlabel)

        return contour_higgs

all_data = load_all("continuum/offset_20/N_k=20000/T=0.0", "resolvents.json.gz").query("k_F == 4.25")

##########################
#####       g        #####
##########################

tasks = [
    (all_data.query("coulomb_scaling == 0 & lambda_screening == 0 & omega_D == 10 & g > 0.25 & g < 3.5"),     "g", r"$g$"),
    (all_data.query("coulomb_scaling == 1 & lambda_screening == 1 & omega_D == 10 & g > 0.7 & g < 3.5"),      "g", r"$g$"),
    (all_data.query("coulomb_scaling == 1 & lambda_screening == 0.0001 & omega_D == 10 & g > 0.7 & g < 3.5"), "g", r"$g$"),
]

fig, axes = plt.subplots(nrows=2, ncols=len(tasks), figsize=(12.8, 6.4), sharex='col', sharey='row')
fig.subplots_adjust(wspace=0.05, hspace=0.1)

import string
for i, axs in enumerate(axes):
    for j, ax in enumerate(axs):
        ax.annotate(
            f"({string.ascii_lowercase[i]}.{j+1})",
            xy=(0, 1), xycoords='axes fraction', xytext=(+0.5, -0.5), textcoords='offset fontsize', 
            verticalalignment='top', fontfamily='serif', bbox=dict(facecolor='1.0', edgecolor='none', pad=3.0))

contour_for_colorbar = None
for i, (data_query, x_column, xlabel) in enumerate(tasks):
    plotter = HeatmapPlotter(data_query, x_column, xlabel=xlabel)
    contour_for_colorbar = plotter.plot(axes[:, i], labels=not bool(i))

cbar = fig.colorbar(contour_for_colorbar, ax=axes[:, -1], orientation='vertical', fraction=0.1, pad=0.05, extend='max')
cbar.set_label(r'$A$ [$\mathrm{meV}^{-1}$]')

filename = os.path.join(BUILD_DIR, f"{os.path.basename(__file__).split('.')[0]}{FILE_ENDING}")
fig.savefig(filename)

##########################
#####     lambda     #####
##########################
all_data = load_all("continuum/offset_10/N_k=20000/T=0.0", "resolvents.json.gz").query("k_F == 4.25")
tasks = [
    (all_data.query("coulomb_scaling == 1 & omega_D == 10 & g == 0.5 & lambda_screening > 1e-2"), "lambda_screening", r"$\lambda$", "log"),
    #(all_data.query("coulomb_scaling == 1 & omega_D == 10 & g == 0.7 & lambda_screening > 1e-2"), "lambda_screening", r"$\lambda$", "log"),
    #(all_data.query("coulomb_scaling == 1 & lambda_screening == 0.0001 & g == 1 & omega_D < 21"),             "omega_D",          f"omega_D_small_screening{FILE_ENDING}", r"$\omega_\mathrm{D} [\mathrm{meV}]$", r"$\lambda = 0.0001$"),
    #(all_data.query("coulomb_scaling == 1 & lambda_screening == 1 & g == 1 & omega_D < 21"),                  "omega_D",          f"omega_D_large_screening{FILE_ENDING}", r"$\omega_\mathrm{D} [\mathrm{meV}]$", r"$\lambda = 1$"),
    #(all_data.query("coulomb_scaling == 0 & lambda_screening == 0 & g == 1 & omega_D < 21"),                  "omega_D",          f"omega_D_no_coulomb{FILE_ENDING}",      r"$\omega_\mathrm{D} [\mathrm{meV}]$", "No Coulomb"),
]

fig, axes = plt.subplots(nrows=2, ncols=len(tasks), figsize=(6.4, 6.4), sharex='col', sharey='row')
fig.subplots_adjust(wspace=0.05, hspace=0.1)

import string
for i, ax in enumerate(axes):
    ax.annotate(
        f"({string.ascii_lowercase[i]})",
        xy=(0, 1), xycoords='axes fraction', xytext=(+0.5, -0.5), textcoords='offset fontsize', 
        verticalalignment='top', fontfamily='serif', bbox=dict(facecolor='1.0', edgecolor='none', pad=3.0))

contour_for_colorbar = None
for i, (data_query, x_column, xlabel, xscale) in enumerate(tasks):
    plotter = HeatmapPlotter(data_query, x_column, xlabel=xlabel, xscale=xscale)
    contour_for_colorbar = plotter.plot(axes[:], labels=not bool(i))

cbar = fig.colorbar(contour_for_colorbar, ax=axes[:], orientation='vertical', fraction=0.1, pad=0.05, extend='max')
cbar.set_label(r'$A$ [$\mathrm{meV}^{-1}$]')

filename = os.path.join(BUILD_DIR, f"{os.path.basename(__file__).split('.')[0]}_screening{FILE_ENDING}")
fig.savefig(filename)
