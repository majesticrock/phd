import matplotlib.pyplot as plt
import lib.resolvent_peak as rp

reversed = True

T = 0.0
U = 3.7
V = 1.0
name = f"T={T}/U={U}/V={V}"
folder = "data/modes/square/dos_900/"
name_suffix = "AFM"
fig, ax = plt.subplots()

peak = rp.Peak(f"{folder}{name}", name_suffix, initial_search_bounds=(1., 4.))
peak.improved_peak_position(x0_offset=1e-2, gradient_epsilon=1e-12)
popt, pcov, w_space, y_data = peak.fit_real_part(range=0.01, begin_offset=1e-10, reversed=True)

ax.text(0.05, 0.35, f"$a={popt[0]:.5f}$", transform = ax.transAxes)
ax.text(0.05, 0.3, f"$b={popt[1]:.5f}$", transform = ax.transAxes)
ax.plot(w_space, y_data, linewidth=1.75*plt.rcParams["lines.linewidth"], label="Data")
ax.plot(w_space, rp.linear_function(w_space, *popt), "k--", label="Fit")

plt.legend()
ax.set_xlabel(r"$\ln((z - z_0) / t)$")
ax.set_ylabel(r"$\ln(\Re G^\mathrm{ret}(z - z_0))$")
fig.tight_layout()

import os
plt.savefig(f"python/build/{os.path.basename(__file__).split('.')[0]}.pdf")
plt.show()