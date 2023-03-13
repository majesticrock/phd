import matplotlib.pyplot as plt
import numpy as np

nameU = "-2.00"
folder = "T0"

UPPER_LIM = 0.03

sub_folders = ["10", "20", "30", "40", "50", "60"]
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']
x_lin = np.linspace(0, UPPER_LIM, 1000)

fig, axs = plt.subplots(2, 1, sharex=True, figsize=(8, 12))

for idx, sub in enumerate(sub_folders):
    diff = np.loadtxt(f"data/{folder}/U_modes/{sub}/{nameU}_special.txt").transpose()

    diff = diff[abs(diff) < UPPER_LIM]
    print(np.min(np.abs(diff)))
    hist, bins = np.histogram(diff, bins=200, range=(0, UPPER_LIM), density=True)
    from scipy.stats import gaussian_kde
    kde = gaussian_kde(diff)
    probability = kde(x_lin)

    axs[0].hist(diff, bins=bins, density=True, label=f"$L={sub}$", color=colors[idx], alpha=0.4)
    axs[1].plot(x_lin, probability, label=f"$L={sub}$", color=colors[idx])


axs[0].set_ylabel("Count")
axs[1].set_ylabel("Probability density (scipy)")
axs[1].set_xlabel("$|\\Delta E|$")

axs[1].legend()
plt.tight_layout()
import os
plt.savefig(f"python/build/{os.path.basename(__file__).split('.')[0]}_U={nameU}.pdf")
plt.show()
