import matplotlib.pyplot as plt
import numpy as np

nameU = "-2.00"
folder = "T0"

UPPER_LIM = 0.02

sub_folders = [12, 18, 20, 24, 30, 36, 40, 46, 50, 60]
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']
x_lin = np.linspace(0, UPPER_LIM, 1000)

fig, axs = plt.subplots(2, 1, sharex=True, figsize=(8, 12))

from scipy.stats import gaussian_kde
max_prob = np.full((2, len(sub_folders)), 0.0)

for idx, sub in enumerate(sub_folders):
    diff = (np.loadtxt(f"data/{folder}/U_modes/{sub}/{nameU}_special.txt").transpose())

    diff = diff[abs(diff) < UPPER_LIM]
    print(np.min(np.abs(diff)))
    hist, bins = np.histogram(diff, bins=200, density=True)#
    
    kde = gaussian_kde(diff)
    probability = kde(x_lin)

    pos = np.argmax(probability)
    max_prob[0][idx] = x_lin[pos]
    max_prob[1][idx] = probability[pos]
    
    axs[0].hist(diff, bins=bins, density=True, label=f"$L={sub}$", color=colors[idx], alpha=0.4)
    axs[1].plot(x_lin, probability, label=f"$L={sub}$", color=colors[idx])

axs[1].plot(max_prob[0], max_prob[1], color="k", linestyle="-", marker="o", label="Maxima")

axs[0].set_ylabel("Adjusted count")
axs[1].set_ylabel("Probability density (scipy)")
axs[1].set_xlabel("$|\\Delta E|$")

#axs[1].set_ylim(1e-1, np.max(max_prob[1])*1.1)
#axs[1].set_yscale("log")
#axs[1].set_xlim(1e-6, UPPER_LIM)
#axs[1].set_xscale("log")

axs[1].legend()
plt.tight_layout()
import os
plt.savefig(f"python/build/{os.path.basename(__file__).split('.')[0]}_U={nameU}.pdf")
#plt.show()


from scipy.optimize import curve_fit
def func(x, a, n):
    return a * x**n

plt.figure()

popt, pcov = curve_fit(func, sub_folders, max_prob[0], p0=(1, -2))
print(popt, pcov)
plt.plot(sub_folders, max_prob[0], color="k", linestyle="-", marker="o", label="Maxima")
x_lin = np.linspace(np.min(sub_folders), np.max(sub_folders))
plt.plot(x_lin, func(x_lin, *popt), color="orange", linestyle="--", label="Fit")

plt.yscale("log")
plt.xscale("log")

plt.ylabel("$\\Delta E (Maximum)$")
plt.xlabel("$L$")
plt.show()