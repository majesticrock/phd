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

screenings = np.log(screenings)
peak_positions = np.log(peak_positions / gaps)

ax.plot(screenings, peak_positions, "o", label=r"$\omega_P$")

from ez_fit import ez_linear_fit
cut = slice(0, 45)
popt, pcov = ez_linear_fit(screenings[cut], peak_positions[cut], ax, x_bounds=(min(screenings), max(screenings)), label=r"Fit $a \ln(\omega) + b$")
ax.text(0.3, 0.9, f"$a = {popt[0]:1.5f} \pm {pcov[0][0]:1.5f}$", transform=ax.transAxes)
ax.text(0.3, 0.8, f"$b = {popt[1]:1.5f} \pm {pcov[1][1]:1.5f}$", transform=ax.transAxes)
ax.set_xlabel(r'$\ln (\lambda)$')
ax.set_ylabel(r'$\ln (\omega / \Delta$)')


ax.legend()

fig.tight_layout()
plt.savefig("python/continuum/build/phase_peaks.pdf")
plt.show()