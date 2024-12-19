import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate

OMEGA_D = 0.01 # eV
K_F = 4.25
LAMBDA = 1e-4
SCREENING = LAMBDA * 0.4107320221286488672 / np.sqrt(K_F)
def coulomb_fock(x):
    if LAMBDA == 0:
        return 0
    k_diff = x - 1
    k_sum = x + 1
    ln_factor = (SCREENING * SCREENING + 1 - x * x) / (2.0 * x)
    log_expression = np.log((SCREENING * SCREENING + k_sum * k_sum) / (SCREENING * SCREENING + k_diff * k_diff))
    return -0.12652559550141668 * K_F * (
        1 + ln_factor * log_expression + SCREENING * (np.arctan(k_diff / SCREENING) - np.arctan(k_sum / SCREENING))
    )

E_F = 0.5 * K_F**2 + coulomb_fock(1) # eV
RHO_F = K_F / (2 * np.pi**2)
G = 0.5 / RHO_F # M^2 / omega_D, yields Delta_max = 2.3 meV
G_TILDE = (G * OMEGA_D / E_F) / (np.pi**2)
omega_tilde = OMEGA_D / E_F

OFFSET = 0.002
EPS = 1e-8
x = np.linspace(1 - OFFSET, 1 + OFFSET, 1000)

def bare_dispersion(x):
    return 0.5 * K_F**2 * ( x**2 - 1 ) 

def phonon_fock(x):
    sqrt_plus = np.sqrt(x*x + omega_tilde)
    sqrt_minus = np.sqrt(x*x - omega_tilde)

    term1 = sqrt_plus  * np.log(  (1 + sqrt_plus)  / np.abs( (1 - sqrt_plus ) + EPS) )
    term2 = sqrt_minus * np.log(  (1 + sqrt_minus) / np.abs( (1 - sqrt_minus) + EPS) )
    return G_TILDE * (term1 - term2)

print("Phononic Fock energy at singularity =", phonon_fock(np.sqrt(1 - omega_tilde)), "eV")
print("Integral on paper =", integrate.quad(phonon_fock, 0.5, 1.))

plt.plot(x, bare_dispersion(x), label=r"$\epsilon_0$")
plt.plot(x, phonon_fock(x), label=r"$\epsilon_\mathrm{Fock}^\mathrm{Ph}$")
plt.plot(x, bare_dispersion(x) + phonon_fock(x), label=r"$\epsilon$")
plt.plot(x, coulomb_fock(x) - coulomb_fock(1), label=r"$\epsilon_\mathrm{Fock}^\mathrm{C} (x, \lambda=10^{-4})$")
#plt.title(r"$\sqrt{x^2 + \tilde{ \omega }}\ln\left|\frac{1+\sqrt{x^2 + \tilde{ \omega }}}{1-\sqrt{x^2 + \tilde{ \omega }}}\right| - \sqrt{x^2 - \tilde{ \omega }}\ln\left|\frac{1+\sqrt{x^2 - \tilde{ \omega }}}{1-\sqrt{x^2 - \tilde{ \omega }}}\right|$")
plt.xlabel(r"$x = k / k_\mathrm{F}$")
plt.ylabel(r"$\epsilon (x) [eV]$")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
