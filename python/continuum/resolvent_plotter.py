import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
from matplotlib import colors
import __path_appender as __ap
__ap.append()
import continued_fraction_pandas as cf

BUILD_DIR = "python/continuum/build/"
FILE_ENDING = ".pdf"

class HeatmapPlotter:
    def __init__(self, data_frame_param, parameter_name, xlabel='Y-axis', zlabel=r'$A$ [$\mathrm{meV}^{-1}$]', title='Spectral functions', xscale="linear", yscale="linear"):
        data_frame = data_frame_param.sort_values(parameter_name).reset_index(drop=True)
        
        self.y_gap = np.linspace(0., 5., 2000)
        self.y_mev = np.linspace(0., 60., 2000)
        self.x = (data_frame[parameter_name]).to_numpy()
        self.resolvents = [cf.ContinuedFraction(pd_row, messages=False) for index, pd_row in data_frame.iterrows()]
        self.gaps = [2e-3 * gap for gap in data_frame["Delta_max"]]

        self.xlabel = xlabel
        self.zlabel = zlabel
        self.title = title
        
        self.xscale = xscale
        self.yscale = yscale
        
        self.plot()

    def plot(self, cmap='viridis'):
        self.fig, self.axes = plt.subplots(2, 2, sharex="col", sharey="row", figsize=(8, 6))
        self.fig.subplots_adjust(wspace=0, hspace=0.1)
        
        spectral_functions_higgs     = np.array([res.spectral_density(gap  * self.y_gap + 2e-5j, "amplitude_SC") for res, gap in zip(self.resolvents, self.gaps)]).transpose()
        spectral_functions_phase     = np.array([res.spectral_density(gap  * self.y_gap + 2e-5j, "phase_SC") for res, gap in zip(self.resolvents, self.gaps)]).transpose()
        spectral_functions_higgs_mev = np.array([res.spectral_density(1e-3 * self.y_mev + 2e-5j, "amplitude_SC") for res in self.resolvents]).transpose()
        spectral_functions_phase_mev = np.array([res.spectral_density(1e-3 * self.y_mev + 2e-5j, "phase_SC") for res in self.resolvents]).transpose()
        
        vmax = max(spectral_functions_higgs.max(), spectral_functions_higgs_mev.max(), spectral_functions_phase.max(), spectral_functions_phase_mev.max())
        
        levels = np.linspace(0., min(1.5, vmax), 51, endpoint=True)
        
        contour = self.axes[0][0].contourf(self.x, self.y_gap, spectral_functions_higgs,     cmap=cmap, levels=levels, extend='max')
        _c1____ = self.axes[0][1].contourf(self.x, self.y_gap, spectral_functions_phase,     cmap=cmap, levels=levels, extend='max')
        _c2____ = self.axes[1][0].contourf(self.x, self.y_mev, spectral_functions_higgs_mev, cmap=cmap, levels=levels, extend='max')
        _c3____ = self.axes[1][1].contourf(self.x, self.y_mev, spectral_functions_phase_mev, cmap=cmap, levels=levels, extend='max')
        
        contour.set_edgecolor('face')
        _c1____.set_edgecolor('face')
        _c2____.set_edgecolor('face')
        _c3____.set_edgecolor('face')
        
        cbar = self.fig.colorbar(contour, ax=self.axes, orientation='vertical', fraction=0.046, pad=0.04, extend='max')
        cbar.set_label(self.zlabel)

        self.axes[0][0].set_ylabel(r"$\omega [2 \Delta]$")
        self.axes[1][0].set_ylabel(r"$\omega [\mathrm{meV}]$")
        
        self.axes[1][0].set_xlabel(self.xlabel)
        self.axes[1][1].set_xlabel(self.xlabel)
        
        self.axes[0][0].set_title(r"$\mathcal{A}_\mathrm{Higgs}$")
        self.axes[0][1].set_title(r"$\mathcal{A}_\mathrm{Phase}$")
        
        self.axes[0][0].set_xscale(self.xscale)
        self.axes[0][1].set_xscale(self.xscale)
        self.axes[0][0].set_yscale(self.yscale)
        self.axes[1][0].set_yscale(self.yscale)
        self.fig.suptitle(self.title)

    def show_plot(self):
        plt.show()

    def save_plot(self, filename):
        self.fig.savefig(filename)


import __path_appender as __ap
__ap.append()
from get_data import *
import concurrent.futures

def load_data():
    return load_all("continuum/offset_20/N_k=20000/T=0.0", "resolvents.json.gz").query("k_F == 4.25")

# Function to create and save a heatmap
def create_and_save_heatmap(data_query, x_column, filename, xlabel, title, xscale="linear"):
    plotter = HeatmapPlotter(data_query, x_column, xlabel=xlabel, title=title, xscale=xscale)
    plotter.save_plot(os.path.join(BUILD_DIR, filename))

if __name__ == "__main__":
    all_data = load_data()

    # Define tasks for heatmap creation and saving
    tasks = [
        (all_data.query("coulomb_scaling == 1 & omega_D == 10 & g == 0.5 & lambda_screening > 1e-2"),             "lambda_screening", f"screening_g05{FILE_ENDING}",           r"$\lambda$",                          r"$g = 0.5$", "log"),
        (all_data.query("coulomb_scaling == 1 & omega_D == 10 & g == 0.7 & lambda_screening > 1e-2"),             "lambda_screening", f"screening_g07{FILE_ENDING}",           r"$\lambda$",                          r"$g = 0.7$", "log"),
        (all_data.query("coulomb_scaling == 1 & lambda_screening == 0.0001 & g == 1 & omega_D < 21"),             "omega_D",          f"omega_D_small_screening{FILE_ENDING}", r"$\omega_\mathrm{D} [\mathrm{meV}]$", r"$\lambda = 0.0001$"),
        (all_data.query("coulomb_scaling == 1 & lambda_screening == 1 & g == 1 & omega_D < 21"),                  "omega_D",          f"omega_D_large_screening{FILE_ENDING}", r"$\omega_\mathrm{D} [\mathrm{meV}]$", r"$\lambda = 1$"),
        (all_data.query("coulomb_scaling == 0 & lambda_screening == 0 & g == 1 & omega_D < 21"),                  "omega_D",          f"omega_D_no_coulomb{FILE_ENDING}",      r"$\omega_\mathrm{D} [\mathrm{meV}]$", "No Coulomb"),
        (all_data.query("coulomb_scaling == 1 & lambda_screening == 0.0001 & omega_D == 10 & g > 0.7 & g < 3.5"), "g",                f"g_small_screening{FILE_ENDING}",       r"$g$",                                 r"$\lambda = 0.0001$"),
        (all_data.query("coulomb_scaling == 1 & lambda_screening == 1 & omega_D == 10 & g > 0.7 & g < 3.5"),      "g",                f"g_large_screening{FILE_ENDING}",       r"$g$",                                 r"$\lambda = 1$"),
        (all_data.query("coulomb_scaling == 0 & lambda_screening == 0 & omega_D == 10 & g > 0.25 & g < 3.5"),     "g",                f"g_no_coulomb{FILE_ENDING}",            r"$g$",                                 "No Coulomb"),
    ]

    # Use ProcessPoolExecutor to parallelize the heatmap generation and saving
    with concurrent.futures.ProcessPoolExecutor() as executor:
        # Map the tasks to executor
        futures = [executor.submit(create_and_save_heatmap, *task) for task in tasks]

        # Ensure all tasks are completed
        concurrent.futures.wait(futures)