import numpy as np
import matplotlib.pyplot as plt
import lib.continued_fraction as cf
from lib.iterate_containers import iterate_containers
from lib.extract_key import *
# Calculates the resolvent in w^2

prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

Ts = np.array([0.])
Us = np.array([-2.])
Vs = np.array(["50.0", "25.0", "15.0", "10.0", "8.0", "6.0", "4.0", "3.5", "3.0", "2.5", "2.0", "1.5", 
               "1.4", "1.3", "1.2", "1.1", "1.0", "0.9", "0.8", "0.7", "0.6", "0.5", 
               "0.45", "0.4", "0.35", "0.3", "0.25", "0.2", "0.15", "0.1", "0.07", 
               "0.05", "0.04", "0.03", "0.02", "0.01",
               "0.007", "0.005", "0.004", "0.003", "0.002", "0.0015", "0.001", 
               "0.0007", "0.0005", "0.0004", "0.0003", "0.0002", "0.00015", "0.0001", 
               "0.00007", "0.00005", "0.00003"
               ])

folder = "data/modes/square/dos_900/"
name_suffix = "phase_SC"

weights = np.zeros(len(Vs))
counter = 0

#fig_check, ax_check = plt.subplots()
for T, U, V in iterate_containers(Ts, Us, Vs):
    name = f"T={T}/U={U}/V={V}"
    lower = float(V)
    upper = 8 * float(V) + 2
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
    scipy_result = opt.fmin_l_bfgs_b(min_func, search_bounds[1] - 2e-1, bounds=[search_bounds], approx_grad=True, epsilon=1e-12, factr=1e7)

    if scipy_result[2]["warnflag"] != 0:
        print(f"We might not have found the peak for V={V}!\nWe found ", peak_pos_value, " and\n", scipy_result)
    peak_pos_value = scipy_result[0][0]
    
    data, data_real, w_lin, res = cf.resolvent_data_log_z(f"{folder}{name}", name_suffix, lower_edge=peak_pos_value,
                                                              range=0.01, begin_offset=1e-8,
                                                              number_of_values=2000, imaginary_offset=1e-6, xp_basis=True, messages=False)

    fit_data = np.log(np.abs(data_real))

    from scipy.optimize import curve_fit
    def func(x, a, b):
        return a * x + b

    popt, pcov = curve_fit(func, w_lin, fit_data)
    #line, = ax_check.plot(w_lin, fit_data, ls="-", label=f"V={V}")
    #ax_check.plot(w_lin, func(w_lin, *popt), ls="--", linewidth=4, color=line.get_color())
    weights[counter] = popt[1]

    counter += 1

v_data = np.log(np.array([float(v) for v in Vs]))
cut = -18
fig, ax = plt.subplots()

ax.plot(v_data[cut:], weights[cut:], "X", label="Fitted data")
ax.plot(v_data[:cut], weights[:cut], "o", label="Omitted data")

popt, pcov = curve_fit(func, v_data[cut:], weights[cut:])
v_lin = np.linspace(v_data.min(), v_data.max(), 500)
ax.plot(v_lin, func(v_lin, *popt), label="Fit")

ax.text(0.05, 0.35, f"$a={popt[0]:.5f}$", transform = ax.transAxes)
ax.text(0.05, 0.3, f"$b={popt[1]:.5f}$", transform = ax.transAxes)

ax.set_xlabel(r"$\ln(V / t)$")
ax.set_ylabel(r"$\ln(w_0 \cdot t)$")
ax.legend()
fig.tight_layout()

import os
plt.savefig(f"python/build/{os.path.basename(__file__).split('.')[0]}.pdf")
plt.show()
