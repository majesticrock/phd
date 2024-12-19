import numpy as np
from scipy import integrate
from uncertainties import ufloat

DELTA = 0.0002
OMEGA = 0.001

def n_sc(x, delta):
    return 0.5 * x**2 * ( 1 - (x**2 - 1) / np.sqrt((x**2 - 1)**2 + delta**2))

INT_N_FS = 1. / 3.
INT_N_SC = ufloat(*integrate.quad(n_sc, 0, np.inf, DELTA))

PREFACTOR = 1 / np.pi**2
kF = 4.25

print("Hartree:", kF**3 * PREFACTOR * (INT_N_FS - INT_N_SC))


def fock(x, pole, delta, fermi_sea):
    #return -x**2 / (x**2 - y**2 + omega) * ( np.heaviside(1 - x, 0.5) - 0.5 * (1 - (x**2 - 1) / np.sqrt((x**2 - 1)**2 + delta**2)) )
    n_expec = np.heaviside(1 - x, 0.5) - (0.5 * (1 - (x**2 - 1) / np.sqrt((x**2 - 1)**2 + delta**2)) if not fermi_sea else 0.0 )
    return (x**2 / (x + pole)) * n_expec

def singularities(y, omega):
    return (np.sqrt(y**2 - omega), np.sqrt(y**2 + omega))

def test(x):
    return 1 
print("test:", integrate.quad(test, -1, 1, weight='cauchy', wvar=0, limit=200))

FOCK_FACTOR = 0.01
def compute_fock_energy(ys, fermi_sea):
    ret_data = []
    for y in ys:
        POLES = singularities(y, OMEGA)
        INT_FOCK1 = ufloat(*integrate.quad(fock, 0, 2, args=(POLES[0], DELTA, fermi_sea), weight='cauchy', wvar=POLES[0], limit=500)) + ufloat(*integrate.quad(fock, 2, 1e6, args=(POLES[0], DELTA, fermi_sea), weight='cauchy', wvar=POLES[0], limit=200))
        INT_FOCK2 = ufloat(*integrate.quad(fock, 0, 2, args=(POLES[1], DELTA, fermi_sea), weight='cauchy', wvar=POLES[1], limit=500)) + ufloat(*integrate.quad(fock, 2, 1e6, args=(POLES[1], DELTA, fermi_sea), weight='cauchy', wvar=POLES[1], limit=200))
        ret_data.append(FOCK_FACTOR * (INT_FOCK1 - INT_FOCK2))
    return ret_data

print("Fock at k_F:", compute_fock_energy([1.], False))

import matplotlib.pyplot as plt

y_lin = np.linspace(0.995, 1.005, 300)

fock_data = compute_fock_energy(y_lin, False)
plot_data = np.array([eps.n for eps in fock_data])
plot_err = 1.96 * np.array([eps.s for eps in fock_data])
plt.plot(y_lin, plot_data, label="SC phase")
plt.fill_between(y_lin, (plot_data - plot_err), (plot_data + plot_err), alpha=0.1)

fock_data = compute_fock_energy(y_lin, True)
plot_data = np.array([eps.n for eps in fock_data])
plot_err = 1.96 * np.array([eps.s for eps in fock_data])
plt.plot(y_lin, plot_data, label="Fermi sea")
plt.fill_between(y_lin, (plot_data - plot_err), (plot_data + plot_err), alpha=0.1)

plt.legend()
plt.xlabel(r"$k / k_\mathrm{F}$")
plt.ylabel(r"$\epsilon_\mathrm{Fock}^\mathrm{Phonon}$")
plt.tight_layout()
plt.show()