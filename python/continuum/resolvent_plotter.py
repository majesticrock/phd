import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
from matplotlib import colors
import __path_appender as __ap
__ap.append()
import continued_fraction_pandas as cf

class HeatmapPlotter:
    def __init__(self, data_frame, parameter_name, xlabel=r'$\omega$ $[2 \Delta]$', ylabel='Y-axis', zlabel=r'$A$ [$\mathrm{meV}^{-1}$]'):
        self.x = np.linspace(0., 5., 2000)
        self.y = (data_frame[parameter_name]).to_numpy()
        self.resolvents = [cf.ContinuedFraction(pd_row, messages=False) for index, pd_row in data_frame.iterrows()]
        self.gaps = [2e-3 * gap for gap in data_frame["Delta_max"]]

        self.xlabel = xlabel
        self.ylabel = ylabel
        self.zlabel = zlabel
        self.fig, self.ax = plt.subplots()
        
        self.plot()

    def plot(self, resolvent_name="amplitude_SC", cmap='viridis'):
        spectral_functions = np.array([res.spectral_density(gap * self.x + 1e-3j, resolvent_name) for res, gap in zip(self.resolvents, self.gaps)])
        
        #y_fine = np.linspace(self.y.min(), self.y.max(), 200)
        #Z_interp = np.zeros((len(y_fine), len(self.x)))
        #for i in range(len(self.x)):
            #interp_func = interp1d(self.y, spectral_functions[:, i], kind='nearest')
            #Z_interp[:, i] = interp_func(y_fine)
        
        #divnorm=colors.TwoSlopeNorm(vmin=spectral_functions.min(), vcenter=0., vmax=spectral_functions.max())
        contour = self.ax.pcolormesh(self.x, self.y, spectral_functions, cmap=cmap)
        
        cbar = self.fig.colorbar(contour, ax=self.ax)
        cbar.set_label(self.zlabel)

        self.ax.set_xlabel(self.xlabel)
        self.ax.set_ylabel(self.ylabel)

    def customize_ticks(self, x_ticks=None, y_ticks=None):
        if x_ticks is not None:
            self.ax.set_xticks(x_ticks)
        if y_ticks is not None:
            self.ax.set_yticks(y_ticks)

    def show_plot(self):
        plt.show()

    def save_plot(self, filename):
        self.fig.savefig(filename)


import __path_appender as __ap
__ap.append()
from get_data import *

all_data = load_all("continuum/offset_10/N_k=20000/T=0.0", "resolvents.json.gz")

#g_small_screening = HeatmapPlotter(all_data.query("coulomb_scaling == 1 & lambda_screening == 0.0001 & k_F == 4.25 & omega_D == 10"), "g", ylabel=r"$g$")
#g_large_screening = HeatmapPlotter(all_data.query("coulomb_scaling == 1 & lambda_screening == 1      & k_F == 4.25 & omega_D == 10"), "g", ylabel=r"$g$")
g_no_coulomb      = HeatmapPlotter(all_data.query("coulomb_scaling == 0 & lambda_screening == 0      & k_F == 4.25 & omega_D == 10"), "g", ylabel=r"$g$")

print(all_data.query("coulomb_scaling == 0 & lambda_screening == 0      & k_F == 4.25 & omega_D == 10")["g"])

plt.show()