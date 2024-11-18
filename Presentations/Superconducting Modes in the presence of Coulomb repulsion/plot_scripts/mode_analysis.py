import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

XTYPE= "g"

df = pd.read_pickle("modes/higgs_0.pkl").sort_values("g")
filtered_df = df[df['energies'].apply(len) > 0].reset_index()
fig, ax = plt.subplots()

for i in range(4):
    first_entries = filtered_df["energies"].apply(lambda x: x[i] if len(x) > i else None)
    valid_idx = first_entries.first_valid_index()  # Get the first index where first_entries is not None

    if valid_idx is not None:  # Ensure there's a valid entry to plot
        ax.plot(filtered_df[XTYPE], first_entries, color=f"C{i}")
        ax.scatter(filtered_df[XTYPE].iloc[valid_idx], first_entries.iloc[valid_idx], color=f"C{i}", marker="*")
        ax.axvline(filtered_df[XTYPE].iloc[valid_idx], ymax=first_entries.iloc[valid_idx], color=f"C{i}", ls="--")

ax.set_xlabel(r"$\Delta_\mathrm{max}$ $[\mathrm{meV}]$")
ax.set_ylabel(r"$\omega$ $[\mathrm{meV}]$")
ax.set_xlim(0, filtered_df[XTYPE].max())

fig.tight_layout()
fig.savefig(f"plots/{os.path.basename(__file__).split('.')[0]}.pdf")