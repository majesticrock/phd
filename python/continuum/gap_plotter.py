import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
from matplotlib import colors

class HeatmapPlotter:
    def __init__(self, data_frame, parameter_name, xlabel=r'$k$ $[k_\mathrm{F}]$', ylabel='Y-axis', zlabel=r'$\Delta$ [meV]'):
        data_size = int(data_frame["inner_discretization"].iloc[0] / 2)
        n_skip = 4
        self.x = (data_frame["data"].iloc[0]["ks"] / data_frame["k_F"].iloc[0]).to_numpy()[data_size:3*data_size:n_skip]
        self.y = (data_frame[parameter_name]).to_numpy()
        self.z = np.array([np.add(row['Delta_Coulomb'], row['Delta_Phonon'])[data_size:3*data_size:n_skip] for row in ( data_frame['data'] )])
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.zlabel = zlabel
        self.fig, self.ax = plt.subplots()
        
        self.plot()

    def plot(self, cmap='seismic'):
        y_fine = np.linspace(self.y.min(), self.y.max(), 200)
        Z_interp = np.zeros((len(y_fine), len(self.x)))
        for i in range(len(self.x)):
            interp_func = interp1d(self.y, self.z[:, i], kind='cubic')
            Z_interp[:, i] = interp_func(y_fine)
        
        divnorm=colors.TwoSlopeNorm(vmin=-self.z.max(), vcenter=0., vmax=self.z.max())
        contour = self.ax.pcolormesh(self.x, y_fine, Z_interp, cmap=cmap, norm=divnorm)
        
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

all_data = load_all("continuum/offset_10/N_k=20000/T=0.0", "gap.json.gz")

g_small_screening = HeatmapPlotter(all_data.query("coulomb_scaling == 1 & lambda_screening == 0.0001 & k_F == 4.25 & omega_D == 10"), "g", ylabel=r"$g$")
g_large_screening = HeatmapPlotter(all_data.query("coulomb_scaling == 1 & lambda_screening == 1      & k_F == 4.25 & omega_D == 10"), "g", ylabel=r"$g$")
g_no_coulomb      = HeatmapPlotter(all_data.query("coulomb_scaling == 0 & lambda_screening == 0      & k_F == 4.25 & omega_D == 10"), "g", ylabel=r"$g$")


plt.show()