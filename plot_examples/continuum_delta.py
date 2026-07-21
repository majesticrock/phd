import matplotlib.pyplot as plt
from mrock.get_data import *

# Create a DataLoader instance.
# By default, this will look for data in the package's default data directory
# or in the directory specified by the MROCK_DATA_DIR environment variable.
data_loader = DataLoader()

# Plot window around k_F and coupling strength used for both datasets.
X_BOUNDS = [1 - 0.004, 1 + 0.004]
G = 0.15

fig, ax = plt.subplots()

# Load the "exact" continuum gap data and normalize k by k_F.
main_df = data_loader.load_panda(
    "continuum",
    "offset_10",
    "gap.json.gz",
    **continuum_params(
        N_k=20000,
        T=0,
        coulomb_scaling=0,
        screening=0,
        k_F=4.25,
        g=G,
        omega_D=10
    )
)
pd_data = main_df["data"]
pd_data["ks"] /= main_df["k_F"]

# Restrict to a range around the plotted X_BOUNDS and plot the exact gap.
plot_range = pd_data.query(f'ks > {0.8 * X_BOUNDS[0]} & ks < {1.2 * X_BOUNDS[1]}')
plot_range.plot("ks", "Delta_Phonon", ax=ax, label=r"$\Delta_\mathrm{exact}$", ls="-", c="C0")

# Load the approximate-theta continuum gap data, normalize k and plot positive gaps.
main_df = data_loader.load_panda(
    "continuum",
    "theta_approx",
    "gap.json.gz",
    **continuum_params(
        N_k=8000,
        T=0,
        coulomb_scaling=0,
        screening=0,
        k_F=4.25,
        g=G,
        omega_D=10
    )
)
pd_data = main_df["data"]
pd_data["ks"] /= main_df["k_F"]

# Plot only entries with positive gap from the approximation.
plot_range = pd_data.query(f'Delta_Phonon > 0')
plot_range.plot("ks", "Delta_Phonon", ax=ax, label=r"$\Delta_\mathrm{approx}$", ls="--", c="C1")

# Mark the approximate gap edges visually with dashed vertical/horizontal connectors.
delta_edges = [plot_range["ks"].iloc[0], plot_range["ks"].iloc[-1]]
ax.plot([X_BOUNDS[0], delta_edges[0], delta_edges[0]],
        [0, 0, plot_range["Delta_Phonon"].iloc[0]],
        color="C1", ls="--")
ax.plot([X_BOUNDS[1], delta_edges[1], delta_edges[1]],
        [0, 0, plot_range["Delta_Phonon"].iloc[0]],
        color="C1", ls="--")

# Finalize axis limits, labels, legend and layout.
ax.set_xlim(*X_BOUNDS)
ax.set_xlabel(r"$k / k_\mathrm{F}$")
ax.set_ylabel(r"$\Delta [\mathrm{meV}]$")

ax.legend()
fig.tight_layout()

plt.show()