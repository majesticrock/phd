import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import json
import gzip
import os

def convert_to_numpy(obj):
    if isinstance(obj, list):
        return np.array(obj)
    if isinstance(obj, dict):
        return {k: convert_to_numpy(v) for k, v in obj.items()}
    return obj

file = os.path.join(os.path.dirname(os.path.abspath(__file__)), "chi.json.gz")
rho_F = 0.6591234790556394

with gzip.open(file, 'rt') as f_open:
    jData = json.load(f_open, object_hook=convert_to_numpy)
    
main_df = pd.json_normalize(jData, max_level=1).iloc[0]


tick_labels = [r'$\Gamma$', 'X', 'M', 'R', r'$\Gamma$']
fig, ax = plt.subplots()

chi_vals = [np.concatenate([
    main_df["chi_Gamma_X"][i],
    main_df["chi_X_M"][i],
    main_df["chi_M_R"][i],
    main_df["chi_R_Gamma"][i],
]) for i in range(len(main_df["chi_Gamma_X"]))]

x = np.arange(len(chi_vals[0]))
for o, omega in enumerate(main_df["omegas"]):
    ax.plot(x, chi_vals[o], label=f"$\\omega={omega}$", c=f"C{o}", ls="-")

ax.axhline(rho_F, ls=":", c="k", label=r"$\rho_F$")
# vertical separators
tick_pos = [0]
for i in range(4):
    tick_pos.append(tick_pos[i] + len(main_df["chi_Gamma_X"][0]))
for p in tick_pos:
    ax.axvline(p, color='k', linestyle=':', alpha=0.6)

# x ticks at symmetry points
ax.set_xticks(tick_pos, tick_labels)
ax.set_ylabel(r"$\Re \chi_0(q,\omega)$")
ax.set_xlabel(r"$qa$")
ax.legend()
fig.tight_layout()
plt.show()