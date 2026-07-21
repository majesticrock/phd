import matplotlib.pyplot as plt
from mrock.get_data import *

# Create a DataLoader instance.
# By default, this will look for data in the package's default data directory
# or in the directory specified by the MROCK_DATA_DIR environment variable.
data_loader = DataLoader()

# Prepare figure and plotting variables.
fig, ax = plt.subplots()
N = 16000
SYSTEM = 'bcc'

# Load lattice cut data (gap vs. energy) for the specified lattice/system.
# lattice_cut_params controls sampling (N) and physical parameters like coupling g, interaction U, Fermi energy and Debye freq.
main_df = data_loader.load_panda(
    "lattice_cut",
    f"./{SYSTEM}",
    "gap.json.gz",
    **lattice_cut_params(
        N=N,
        g=1.89,
        U=0.01,
        E_F=-0.5,
        omega_D=0.02
    )
)

# Extract energy axis and gap, then plot.
energy_space = main_df['energies']
delta = main_df['Delta']
ax.plot(energy_space, delta, 'k-')

# Axis labels and layout.
ax.set_xlabel(r'$\epsilon - \mu$')
ax.set_ylabel(r'$\Delta$')

fig.tight_layout()
plt.show()