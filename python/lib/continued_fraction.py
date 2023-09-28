import numpy as np
import gzip
import matplotlib.pyplot as plt

class ContinuedFraction:
    # In Python there is no need to declare one's variables beforehand.
    # How foolish of me to assume otherwise
    
    def terminator(self, w):
        p = w - self.a_infinity
        q = 4 * self.b_infinity**2
        root = np.sqrt(np.real(p**2 - q), dtype=complex)
        return_arr = np.zeros(len(w), dtype=complex)
        for i in range(0, len(w)):
            if(w[i] > self.roots[0]):
                return_arr[i] = (p[i] - root[i]) / (2. * self.b_infinity**2)
            else:
                return_arr[i] = (p[i] + root[i]) / (2. * self.b_infinity**2)
        return return_arr

    def continued_fraction(self, w, withTerminator = True):
        if withTerminator:
            for i in range(0, len(w)):
                if(w[i].real > self.roots[0] and w[i].real < self.roots[1]):
                    w[i] = w[i].real
            G = w - self.A[len(self.A) - self.terminate_at] - self.B[len(self.B) - self.terminate_at] * self.terminator(w)
        else:
            G = w - self.A[len(self.A) - self.terminate_at]
        for j in range(len(self.A) - self.terminate_at - 1, -1, -1):
            G = w - self.A[j] - self.B[j + 1] / G

        return self.B[0] / G
    
    def mark_continuum(self, axes=None):
        if axes is None:
            if self.z_squared:
                plt.axvspan(np.sqrt(self.roots[0]), np.sqrt(self.roots[1]), alpha=.2, color="purple", label="Continuum")
            else:
                plt.axvspan(self.roots[0], self.roots[1], alpha=.2, color="purple", label="Continuum")
        else:
            if self.z_squared:
                axes.axvspan(np.sqrt(self.roots[0]), np.sqrt(self.roots[1]), alpha=.2, color="purple", label="Continuum")
            else:
                axes.axvspan(self.roots[0], self.roots[1], alpha=.2, color="purple", label="Continuum")
    
    def __init__(self, data_folder, filename, z_squared=True):
        self.z_squared = z_squared
        
        file = f"{data_folder}/one_particle.dat.gz"
        with gzip.open(file, 'rt') as f_open:
            one_particle = np.abs(np.loadtxt(f_open).flatten())
        if z_squared:
            self.roots = np.array([np.min(one_particle) * 2, np.max(one_particle) * 2])**2
        else:
            self.roots = np.array([np.min(one_particle) * 2, np.max(one_particle) * 2])
        self.a_infinity = (self.roots[0] + self.roots[1]) * 0.5
        self.b_infinity = (self.roots[1] - self.roots[0]) * 0.25
        
        file = f"{data_folder}/{filename}.dat.gz"
        with gzip.open(file, 'rt') as f_open:
            M = np.loadtxt(f_open)
            self.A = M[0]
            self.B = M[1]
            
        deviation_from_infinity = np.zeros(len(self.A) - 1)
        for i in range(0, len(self.A) - 1):
            deviation_from_infinity[i] = abs((self.A[i] - self.a_infinity) / self.a_infinity) + abs((np.sqrt(self.B[i + 1]) - self.b_infinity) / self.b_infinity)
        self.terminate_at = len(self.A) - np.argmin(deviation_from_infinity) - 1
        print("Terminating at i=", np.argmin(deviation_from_infinity))