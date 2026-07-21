import matplotlib.pyplot as plt
import numpy as np

from scipy.optimize import curve_fit

from mrock.get_data import *
# Create a DataLoader instance.
# By default, this will look for data in the package's default data directory
# or in the directory specified by the MROCK_DATA_DIR environment variable.
data_loader = DataLoader()

# System and sampling parameters for the T_c calculation.
SYSTEM = 'bcc'
N = 10000
params = lattice_cut_params(
    N=N,
    g=2.5,
    U=0.0,
    E_F=-0.5,
    omega_D=0.02
)

# Load the temperature-dependent maximum gap data.
main_df = data_loader.load_panda("lattice_cut", f"./T_C/{SYSTEM}", "T_C.json.gz", **params)

fig, ax = plt.subplots()

# Linear model for fitting Δ^2 vs T near T_c: Δ^2(T) = m T + b.
def linear_model(T, m, b):
    return m*T + b

# Extract temperatures and maximum gaps from the dataset.
Ts = main_df['temperatures']
deltas = np.array(main_df['max_gaps'])

# Choose a fitting window near the transition (exclude low-T points).
cut = np.min([np.argmin(np.abs(Ts - 0.95 * Ts[-1])), len(Ts) - 5])
T_fit = Ts[cut:]
y_fit = (deltas[cut:])**2

# Fit the linearized data to obtain slope m and intercept b.
popt, pcov = curve_fit(linear_model, T_fit, y_fit)
m, b = popt
dm, db = np.sqrt(np.diag(pcov))

# Convert fit parameters to A and Tc using Δ^2 = -A^2 (T - Tc) => A = sqrt(-m), Tc = b/(-m).
A = np.sqrt(-m)
Tc = b / (-m)

# Propagate errors for A and Tc from the covariance matrix.
dA = abs(0.5/np.sqrt(-m) * dm)
dTc = np.sqrt((b/m**2)**2 * dm**2 + (1/m**2)*db**2 + 2*(b/m**2)*(-1/m)*pcov[0,1])

print(f"\nA = {A} ± {dA}")
print(f"Tc = {Tc} ± {dTc}")

# Plot data and the square-root of the linear fit for visual check.
t_lin = np.linspace(T_fit.min(), Tc, 300)
plt.plot(Ts, deltas, "x--", label="data")
plt.plot(t_lin, np.sqrt(linear_model(t_lin, *popt)), "--", label="linearized fit")
plt.axvline(T_fit[0], linestyle=":", color="k")

ax.legend()
ax.set_xlabel(r'$T$')
ax.set_ylabel(r'$\Delta_\mathrm{max}(T)$')
fig.tight_layout()
plt.show()