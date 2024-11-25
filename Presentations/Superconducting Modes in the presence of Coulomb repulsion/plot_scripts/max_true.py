import matplotlib.pyplot as plt
import os
from __all_data_pickler import load_pickled
import numpy as np

G_MAX_PLOT = 3.4
PICK_G = 2.175

fig, axes = plt.subplots(ncols=3, sharey=True)

all_data = load_pickled()
queries = [
    all_data.query(f"coulomb_scaling == 0 & lambda_screening == 0      & omega_D == 10 & g >= 1 & g <= {G_MAX_PLOT} & Delta_max > 0").sort_values("g"),
    all_data.query(f"coulomb_scaling == 1 & lambda_screening == 1      & omega_D == 10 & g >= 1 & g <= {G_MAX_PLOT} & Delta_max > 0").sort_values("g"),
    all_data.query(f"coulomb_scaling == 1 & lambda_screening == 0.0001 & omega_D == 10 & g >= 1 & g <= {G_MAX_PLOT} & Delta_max > 0").sort_values("g"),
]

for ax, query in zip(axes, queries):
    y_data = np.array([0.5e3 * boundaries[0] for boundaries in query["continuum_boundaries"]])
    ax.plot(query["g"], y_data / query["Delta_max"], "x")
    ax.set_xlabel(r"$g$")
    ax.grid(True)
    if PICK_G is not None:
        ax.scatter(PICK_G, query.query(f"g=={PICK_G}")["continuum_boundaries"].iloc[0][0] * 0.5e3 / query.query(f"g=={PICK_G}")["Delta_max"].iloc[0], color="red")

axes[0].set_ylabel(r"$\Delta_\mathrm{true} / \Delta_\mathrm{max}$")
fig.subplots_adjust(wspace=0.2, hspace=0.1)
fig.savefig(f"plots/{os.path.basename(__file__).split('.')[0]}.pdf")
