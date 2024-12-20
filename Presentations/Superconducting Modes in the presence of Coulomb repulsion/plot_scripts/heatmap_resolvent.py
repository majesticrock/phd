import matplotlib.pyplot as plt
import numpy as np
import __path_appender as __ap
__ap.append()
import continued_fraction_pandas as cf
import os
from get_data import *
from scipy.signal import find_peaks
from __all_data_pickler import DATA_CUTS, load_pickled

BUILD_DIR = "plots/"
FILE_ENDING = ".pdf"
G_MAX_LOAD = 4.25
G_MAX_PLOT = 3.4

class HeatmapPlotter:
    def __init__(self, data_frame_param, parameter_name, xlabel='Y-axis', zlabel=r'$A$ [$\mathrm{eV}^{-1}$]', xscale="linear", yscale="linear"):
        self.data_frame = data_frame_param.sort_values(parameter_name).reset_index(drop=True)
        
        self.y = np.linspace(0., 55., 10000)
        self.x = (self.data_frame[parameter_name]).to_numpy()
        self.resolvents = [cf.ContinuedFraction(pd_row, messages=False, ignore_first=30, ignore_last=88) for index, pd_row in self.data_frame.iterrows()]
        self.max_gaps = [2 * gap for gap in self.data_frame["Delta_max"]]
        self.true_gaps = [1e3 * float(t_gap[0]) for t_gap in self.data_frame["continuum_boundaries"]]
        self.N_data = len(self.max_gaps)
        
        self.g_cuts = np.zeros(len(DATA_CUTS))
        for i in range(len(DATA_CUTS)):
            filtered_df = self.data_frame[self.data_frame['Delta_max'] < DATA_CUTS[i]]
            if len(filtered_df) == 0:
                self.g_cuts[i] = 0
            else:
                closest_row = filtered_df.loc[(DATA_CUTS[i] - filtered_df['Delta_max']).idxmin()]
                self.g_cuts[i] = closest_row['g']

        self.xlabel = xlabel
        self.zlabel = zlabel
        self.xscale = xscale
        self.yscale = yscale

    def identify_modes(self, spectral, pos):
        if self.max_gaps[pos] == 0:
            return np.array([])
        indizes = find_peaks(spectral, distance=int(1 / (self.y[1] - self.y[0])))[0]
        positions = np.array([self.y[i] for i in indizes])
        positions = positions[positions < self.true_gaps[pos]]
        return positions

    def __get_error__(self, key, i):
        __OFFSET__ = 2
        if key == "true_gap":
            return [self.true_gaps[i]     - self.true_gaps[i - __OFFSET__] if i - __OFFSET__ >= 0                   else None,
                    self.true_gaps[i + __OFFSET__] - self.true_gaps[i]     if i + __OFFSET__ <  len(self.true_gaps) else None]
    
        return [self.data_frame[key].iloc[i]     - self.data_frame[key].iloc[i - __OFFSET__] if i - __OFFSET__ >= 0                         else None,
                self.data_frame[key].iloc[i + __OFFSET__] - self.data_frame[key].iloc[i]     if i + __OFFSET__ <  len(self.data_frame["g"]) else None]

    def plot(self, axes, cmap='inferno', labels=True):
        spectral_functions_higgs = np.array([res.spectral_density(1e-3 * self.y + 1e-6j, "amplitude_SC") for res in self.resolvents]).transpose()
        spectral_functions_phase = np.array([res.spectral_density(1e-3 * self.y + 1e-6j, "phase_SC") for res in self.resolvents]).transpose()

        self.HiggsModes = pd.DataFrame([ {
                "resolvent_type": "Higgs",
                "energies": self.identify_modes(spectral_functions_higgs[:, i], i),
                "Delta_max": self.data_frame["Delta_max"].iloc[i],
                "true_gap": self.true_gaps[i],
                "g": self.data_frame["g"].iloc[i],
                "error_g": self.__get_error__("g", i),
                "error_Delta_max": self.__get_error__("Delta_max", i),
                "error_true_gap": self.__get_error__("true_gap", i),
                "T": self.data_frame["T"].iloc[i],
                "omega_D": self.data_frame["omega_D"].iloc[i],
                "E_F": self.data_frame["E_F"].iloc[i],
                "k_F": self.data_frame["k_F"].iloc[i],
                "lambda_screening": self.data_frame["lambda_screening"].iloc[i]
            } for i in range(self.N_data) ])
        self.PhaseModes = pd.DataFrame([ {
                "resolvent_type": "Phase",
                "energies": self.identify_modes(spectral_functions_phase[:, i], i),
                "Delta_max": self.data_frame["Delta_max"].iloc[i],
                "true_gap": self.true_gaps[i],
                "g": self.data_frame["g"].iloc[i],
                "error_g": self.__get_error__("g", i),
                "error_Delta_max": self.__get_error__("Delta_max", i),
                "error_true_gap": self.__get_error__("true_gap", i),
                "T": self.data_frame["T"].iloc[i],
                "omega_D": self.data_frame["omega_D"].iloc[i],
                "E_F": self.data_frame["E_F"].iloc[i],
                "k_F": self.data_frame["k_F"].iloc[i],
                "lambda_screening": self.data_frame["lambda_screening"].iloc[i]
            } for i in range(self.N_data) ])

        vmax = max(spectral_functions_higgs.max(), spectral_functions_phase.max())
        levels = np.linspace(0., min(1.5, vmax), 101, endpoint=True)

        contour_higgs = axes[0].contourf(self.x, self.y, spectral_functions_higgs, cmap=cmap, levels=levels, extend='max', zorder=-20)
        contour_phase = axes[1].contourf(self.x, self.y, spectral_functions_phase, cmap=cmap, levels=levels, extend='max', zorder=-20)
        contour_higgs.set_edgecolor('face')
        contour_phase.set_edgecolor('face')
        
        for ax in axes:
            ax.plot(self.x, self.true_gaps, color="lime", ls=":")
            ax.set_rasterization_zorder(-10)
            ax.set_ylim(0., max(self.y))
            ax.set_xscale(self.xscale)
            ax.set_yscale(self.yscale)

        if(labels):
            axes[0].set_ylabel(r"$\omega [\mathrm{meV}]$")
            axes[1].set_ylabel(r"$\omega [\mathrm{meV}]$")
        axes[1].set_xlabel(self.xlabel)

        return contour_higgs

all_data = load_pickled()

##########################
#####       g        #####
##########################

tasks = [
    (all_data.query(f"coulomb_scaling == 0 & lambda_screening == 0      & omega_D == 10 & g <= {G_MAX_LOAD}"), "g", r"$g$"),
    (all_data.query(f"coulomb_scaling == 1 & lambda_screening == 1      & omega_D == 10 & g <= {G_MAX_LOAD}"), "g", r"$g$"),
    (all_data.query(f"coulomb_scaling == 1 & lambda_screening == 0.0001 & omega_D == 10 & g <= {G_MAX_LOAD}"), "g", r"$g$"),
]

fig, axes = plt.subplots(nrows=2, ncols=len(tasks), figsize=(12.8, 6.4), sharex=True, sharey=True)
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
    
    plotter.HiggsModes.to_pickle(f"modes/higgs_{i}.pkl")
    plotter.PhaseModes.to_pickle(f"modes/phase_{i}.pkl")
    #for j in range(2):
    #    for k in range(len(DATA_CUTS)):
    #        axes[j][i].axvline(plotter.g_cuts[k], color="C4")

cbar = fig.colorbar(contour_for_colorbar, ax=axes[:, -1], orientation='vertical', fraction=0.1, pad=0.05, extend='max')
cbar.set_label(r'$A$ [$\mathrm{meV}^{-1}$]')

axes[0][0].set_xticks([0, 1, 2, 3])
axes[0][0].set_xlim(0, G_MAX_PLOT)
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
