import matplotlib.pyplot as plt
import numpy as np

import numpy as np
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d.proj3d import proj_transform

import matplotlib.colors as colors

class Arrow3D(FancyArrowPatch):
    def __init__(self, x, y, z, dx, dy, dz, *args, **kwargs):
        super().__init__((0, 0), (0, 0), *args, **kwargs)
        self._xyz = (x, y, z)
        self._dxdydz = (dx, dy, dz)

    def draw(self, renderer):
        x1, y1, z1 = self._xyz
        dx, dy, dz = self._dxdydz
        x2, y2, z2 = (x1 + dx, y1 + dy, z1 + dz)

        xs, ys, zs = proj_transform((x1, x2), (y1, y2), (z1, z2), self.axes.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
        super().draw(renderer)
        
    def do_3d_projection(self, renderer=None):
        x1, y1, z1 = self._xyz
        dx, dy, dz = self._dxdydz
        x2, y2, z2 = (x1 + dx, y1 + dy, z1 + dz)

        xs, ys, zs = proj_transform((x1, x2), (y1, y2), (z1, z2), self.axes.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
        return np.min(zs)

def _arrow3D(ax, x, y, z, dx, dy, dz, *args, **kwargs):
    arrow = Arrow3D(x, y, z, dx, dy, dz, *args, **kwargs)
    ax.add_artist(arrow)

def _element(x, i):
    if hasattr(x, "__len__"):
        return x[i]
    return x

# Express the mesh in the cartesian system.
def cartesian(radius, phi):
    return radius * np.cos(phi), radius * np.sin(phi)

def direction(radius_start, phi_start, radius_end, phi_end):
    return (radius_start * np.cos(phi_start) - radius_end * np.cos(phi_end), radius_start * np.sin(phi_start) - radius_end * np.sin(phi_end))

def arrow_terminated_line(ax, r, phi, z, *args, **kwargs):
    ax.plot(*cartesian(r, phi), z, *args, **kwargs)
    kwargs.pop("label")
    _arrow3D(ax, *cartesian(_element(r, 0), _element(phi, 0)), _element(z, 0), 
             *direction(_element(r, 0), _element(phi, 0), _element(r, 1), _element(phi, 1)), 
            _element(z, 0) - _element(z, 1), 
            mutation_scale=20, *args, **kwargs)
    _arrow3D(ax, *cartesian(_element(r, -1), _element(phi, -1)), _element(z, -1), 
             *direction(_element(r, -1), _element(phi, -1), _element(r, -2), _element(phi, -2)), 
            _element(z, - 1) - _element(z, -2), 
            mutation_scale=20, *args, **kwargs)

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(projection='3d')

# Create the mesh in polar coordinates and compute corresponding Z.
r = np.linspace(0, 1.25, 51)
p = np.linspace(0, 2*np.pi, 51)
R, P = np.meshgrid(r, p)

def potential(radius, u, v):
    return 0.5 * u * radius**4 + v * radius**2 - (0.5 if v < 0 else 1.5) * (u**2 / v)

def mininum(u, v):
    return np.sqrt(np.abs(u/v))

Z = potential(R, 1, -1)
X, Y = cartesian(R, P)

# Plot the surface.
def truncate_colormap(cmap, minval=0.0, maxval=0.8, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

ax.plot_surface(X, Y, Z, cmap=truncate_colormap(plt.cm.gist_gray), rstride=2, cstride=2)

# plot modes
phi_length = -4*np.pi/5
arrow_terminated_line(ax, mininum(1, -1), np.linspace(0, phi_length), 0.0, zorder=100, color="C9", label="Phase", linewidth=5)

r_space = np.linspace(mininum(1, -1) - 0.33, mininum(1, -1) + 0.25)
arrow_terminated_line(ax, r_space, phi_length / 2, potential(r_space, 1, -1), zorder=100, color="C1", label="Amplitude", linewidth=5)


# Tweak the limits and add latex math labels.
ax.set_zlim(0, 0.5)
ax.set_xlabel(r'$\mathrm{Re}[ \Delta ]$')
ax.set_ylabel(r'$\mathrm{Im}[ \Delta ]$')
ax.set_zlabel(r'$V(\Delta)$')
ax.legend(loc="center right")

ax.set_xticks([-1, 0, 1])
ax.set_yticks([-1, 0, 1])
ax.set_zticks([0, 0.25, 0.5])
plt.minorticks_off()

fig.tight_layout()
plt.show()