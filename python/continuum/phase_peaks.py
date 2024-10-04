import numpy as np
import matplotlib.pyplot as plt
import __path_appender as __ap
__ap.append()
import continued_fraction_pandas as cf

from get_data import *

main_df = load_all("continuum/offset_10/N_k=8000/T=0.0/coulomb_scaling=1.0", "resolvents.json.gz", condition=["g=0.5"])
print(main_df)
main_df.sort_values("lambda_screening", inplace=True)
main_df.reset_index(inplace=True)


screenings = np.zeros(len(main_df))
peak_positions = np.zeros(len(main_df))

for index, pd_row in main_df.iterrows():
    resolvents = cf.ContinuedFraction(pd_row, messages=False)
    w_lin = np.linspace(-0.005 * pd_row["continuum_boundaries"][1], 1.1 * pd_row["continuum_boundaries"][1], 5000, dtype=complex)
    w_lin += 1e-5j

    y_data = resolvents.spectral_density(w_lin, "phase_SC")
    peak_positions[index] = 1e3 * w_lin.real[np.argmax(y_data)]
    screenings[index] = pd_row['lambda_screening']

fig, ax = plt.subplots()

ax.plot(screenings, peak_positions)

ax.set_xlabel('$\\lambda$')
ax.set_ylabel('$\\omega_P$ [meV]')

fig.tight_layout()
plt.show()