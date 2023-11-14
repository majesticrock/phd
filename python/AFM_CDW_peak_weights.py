import numpy as np
import matplotlib.pyplot as plt
import lib.continued_fraction as cf
from lib.iterate_containers import iterate_containers
from lib.extract_key import *
# Calculates the resolvent in w^2

prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

Ts = np.array([0.])
Us = np.array([3.7, 3.75, 3.8, 3.85, 3.9, 3.925, 3.94, 3.95, 3.97, 3.975, #10 values
               3.985, 3.99, 3.995, 3.997, 3.9985, 3.999, 3.9995, #7
               4.0,
               4.0005, 4.001, 4.0015, 4.003, 4.005, 4.01, 4.015, 4.025, 4.03, 4.05, #10 values
               4.06, 4.075, 4.1, 4.15, 4.2, 4.25, 4.3 #7
               ])
Vs = np.array([1.])

folder = "data/modes/square/dos_900/"
colors = ["orange", "purple"]

name_suffices = ["AFM", "CDW"]
fig, ax = plt.subplots()
u_data = np.array([float(u) for u in Us])

for i, name_suffix in enumerate(name_suffices):
    weights = np.zeros(len(Us))
    counter = 0
    
    for T, U, V in iterate_containers(Ts, Us, Vs):
        name = f"T={T}/U={U}/V={V}"
        
        lower = 1.
        upper = 4.
        data, data_real, w_lin, res = cf.resolvent_data(f"{folder}{name}", name_suffix, lower_edge=lower, upper_edge=upper, 
                                                    number_of_values=20000, imaginary_offset=1e-6, xp_basis=True, messages=False)
        peak_pos_value = w_lin[np.argmax(data)]
        import scipy.optimize as opt
        def min_func(x):
            return res.continued_fraction(x + 1e-6j, True).imag
        
        offset_peak = 0.2
        search_bounds = (0 if peak_pos_value - offset_peak < 0 else peak_pos_value - offset_peak, 
                         np.sqrt(res.roots[0]) if peak_pos_value + offset_peak > np.sqrt(res.roots[0]) else peak_pos_value + offset_peak)
        # epsilon in this function call is *not* the precision but the step size for evaluating the gradient
        # factr is the precision parameter; machine epsilon * factr is the precision we get
        scipy_result = opt.fmin_l_bfgs_b(min_func, search_bounds[1] - 5e-2, bounds=[search_bounds], approx_grad=True, epsilon=1e-12)
        
        if scipy_result[2]["warnflag"] != 0:
            print(f"We might not have found the peak for U={U}!\nWe found ", peak_pos_value, " and\n", scipy_result)
        peak_pos_value = scipy_result[0][0]

        data, data_real, w_lin, res = cf.resolvent_data_log_z(f"{folder}{name}", name_suffix, lower_edge=peak_pos_value,
                                                            range=0.1, begin_offset=1e-10,
                                                            number_of_values=2000, imaginary_offset=1e-6, xp_basis=True, reversed=reversed, messages=False)
        fit_data = np.log(np.abs(data_real))
        from scipy.optimize import curve_fit
        def func(x,  b):
            return - x + b
        popt, pcov = curve_fit(func, w_lin, fit_data)
        
        weights[counter] = np.exp(popt[0])
        counter += 1

    ax.plot(u_data, (weights), marker="X", ls="-", label=f"{name_suffix}", color=colors[i])

ax.set_xlabel(r"$U / t$")
ax.set_ylabel(r"$\ln(w_0 \cdot t)$")
ax.legend()
fig.tight_layout()

import os
plt.savefig(f"python/build/{os.path.basename(__file__).split('.')[0]}.pdf")
plt.show()
