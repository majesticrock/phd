import numpy as np
import matplotlib.pyplot as plt
import __path_appender as __ap
__ap.append()
import continued_fraction_pandas as cf

from get_data import *

main_df = load_all("continuum/offset_10/N_k=20000/T=0.0/coulomb_scaling=1.0", "resolvents.json.gz").query('g == 0.5 & k_F == 4.25 & omega_D == 10')
main_df.sort_values("lambda_screening", inplace=True)
main_df.reset_index(inplace=True)

screenings = np.zeros(len(main_df))
peak_positions = np.zeros(len(main_df))
gaps = np.zeros(len(main_df))

for index, pd_row in main_df.iterrows():
    resolvents = cf.ContinuedFraction(pd_row, messages=False)
    w_lin = np.linspace(-0.005 * pd_row["continuum_boundaries"][1], 
                        1.1 * pd_row["continuum_boundaries"][1] if pd_row["lambda_screening"] > 0.01 else 23., 
                        20000, dtype=complex)
    w_lin += 1e-5j

    y_data = resolvents.spectral_density(w_lin, "phase_SC")
    peak_positions[index] = 1e3 * w_lin.real[np.argmax(y_data)]# / pd_row["Delta_max"]
    screenings[index] = pd_row['lambda_screening']
    gaps[index] = 2 * pd_row["Delta_max"]

fig, ax = plt.subplots()

ax.plot(screenings, peak_positions, "o", label=r"$\omega_P$")
ax.plot(screenings, gaps, "x", label=r"$2 \Delta$")

ax.set_xlabel('$\\lambda$')
ax.set_ylabel('$\\omega$ [meV]')

ax.set_xscale("log")
ax.set_yscale("log")

ax.legend()

fig.tight_layout()
plt.show()