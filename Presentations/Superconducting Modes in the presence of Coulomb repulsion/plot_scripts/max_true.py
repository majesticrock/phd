import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import os
from __all_data_pickler import load_pickled
from get_data import *
import numpy as np
import string

G_MAX_PLOT = 3.3
PICK_G = 2.175

#fig, axes = plt.subplots(2, 3, figsize=(10, 6), gridspec_kw={'width_ratios': [1, 2]}, sharey="row")
fig = plt.figure(figsize=(6.4, 8))
gs = GridSpec(2, 3, hspace=0.3, wspace=0.2, height_ratios=[1.5, 1])
axes_energy = [ fig.add_subplot(gs[0, i]) for i in range(3) ]
ax_gap = fig.add_subplot(gs[1, :])

import matplotlib.colors as mcolors
def adjust_brightness(color, factor):
    return tuple(min(max(c * factor, 0), 1) for c in color)
def rgb_to_grayscale(color):
    return (0.299 * color[0] + 0.587 * color[1] + 0.114 * color[2],) * 3

##### Energy plot

def energy(xi, delta):
    return np.sqrt(xi**2 + 1e-6 * delta**2)

GS = [1., 1.5, 2., 2.5]
SCREENINGS = [0, 1, 1e-4]
OFFS = [0.0014, 0.001, 0.0005]
COLOR_SHIFT = 0.4
linestyles = ["-", "--", "-.", (0, (3, 1, 1, 1, 1, 1))]

last_color = [0, 0, 0]

for i in range(3):
    base_color_rgb = mcolors.to_rgb(plt.cm.tab10(i))
    colors_close_to_C = [
        adjust_brightness(base_color_rgb, 1 + (fac - 1.5) * COLOR_SHIFT) for fac in range(4)
    ]
    last_color[i] = colors_close_to_C[-2]
    for j, g in enumerate(GS):
        main_df = load_panda("continuum", "offset_20", "gap.json.gz", print_date=False,
                            **continuum_params(N_k=20000, T=0, coulomb_scaling=float(i != 0), screening=SCREENINGS[i], k_F=4.25, g=g, omega_D=10))
        pd_data = main_df["data"]
        pd_data["ks"] /= main_df["k_F"]
        if pd_data["Delta_Coulomb"][0] > 0:
            pd_data["Delta_Phonon"] *= -1
            pd_data["Delta_Coulomb"] *= -1
        plot_data = pd_data.query(f"ks > {1-OFFS[i]} & ks < {1+OFFS[i]}")
        axes_energy[i].plot(plot_data["ks"] - 1, 
                            (1e3 * energy(plot_data["xis"] + 1e-3 * plot_data["Delta_Fock"], 
                                         plot_data["Delta_Phonon"] + plot_data["Delta_Coulomb"]) 
                            - main_df["Delta_max"]), 
                            color=colors_close_to_C[j])
        axes_energy[i].set_xlabel(r"$k / k_\mathrm{F} - 1$")
    axes_energy[i].annotate(
            f"({string.ascii_lowercase[i]})",
            xy=(0, 1), xycoords='axes fraction', xytext=(+0.5, -0.5), textcoords='offset fontsize', 
            verticalalignment='top', fontfamily='serif', bbox=dict(facecolor='1.0', edgecolor='black', pad=3))

grayscale_colors = [rgb_to_grayscale(c) for c in colors_close_to_C]
legend_labels = [f"$g={g}$" for g in GS]
legend_handles = [
    plt.Line2D([0], [0], color=grayscale_colors[i], label=legend_labels[i])
    for i in range(len(grayscale_colors))
]

axes_energy[0].set_xticks([-0.001, 0.001])   
axes_energy[1].set_xticks([-0.0007, 0.0007])      
axes_energy[2].set_xticks([-0.0003, 0.0003])    
   
axes_energy[1].set_yticklabels([])      
axes_energy[2].set_yticklabels([])    

axes_energy[0].set_ylabel(r"$E (k) - \Delta_\mathrm{max} [\mathrm{meV}]$")
axes_energy[1].legend(handles=legend_handles, loc="upper center", ncols=2, shadow=True, bbox_to_anchor=(0.5, 1.31), columnspacing=3)

for ax in axes_energy:
    ax.set_ylim(-.99, 1.1)
    

##### Delta_true plot

all_data = load_pickled()
queries = [
    all_data.query(f"coulomb_scaling == 0 & lambda_screening == 0      & omega_D == 10 & g >= 1 & g <= {G_MAX_PLOT} & Delta_max > 0").sort_values("g"),
    all_data.query(f"coulomb_scaling == 1 & lambda_screening == 1      & omega_D == 10 & g >= 1 & g <= {G_MAX_PLOT} & Delta_max > 0").sort_values("g"),
    all_data.query(f"coulomb_scaling == 1 & lambda_screening == 0.0001 & omega_D == 10 & g >= 1 & g <= {G_MAX_PLOT} & Delta_max > 0").sort_values("g"),
]
labels = [
    r"(a)",
    r"(b)",
    r"(c)"
]

for i, (query, label) in enumerate(zip(queries, labels)):
    y_data = np.array([0.5e3 * boundaries[0] for boundaries in query["continuum_boundaries"]])
    ax_gap.plot(query["g"], y_data / query["Delta_max"] - 1, label=label, color=last_color[i], ls=linestyles[i])

ax_gap.legend()
ax_gap.set_xlabel(r"$g$")
ax_gap.set_ylabel(r"$\Delta_\mathrm{true} / \Delta_\mathrm{max} - 1$")
ax_gap.text(3.1, -0.005, "(d)")

#fig.subplots_adjust(wspace=0.2, hspace=0.1)
fig.savefig(f"plots/{os.path.basename(__file__).split('.')[0]}.pdf")
