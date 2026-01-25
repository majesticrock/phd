import numpy as np
import pyvista as pv
import imageio

# Lattice setup
N_X = 24
N_Y = 6

Y_POS_ELECTRON = 3.5
Y_POS_ELECTRON_2 = 4.2
ELECTRON_2_DELAY = 2.0  # seconds after first electron enters

#Video parameters
video_length = 10.0   # seconds
fps = 20              # frames per second
n_frames = int(video_length * fps)
dt = 1 / fps           # time step per frame

electron_velocity = 3.5
coupling_constant = 0.02
scale_factor = 1.5
relaxation_force = 0.03
regularization = 0.1

X, Y = np.meshgrid(np.linspace(0, N_X, N_X, endpoint=False), np.linspace(0., N_Y, N_Y, endpoint=False))
EQUILIBRIUM_LATTICE_POINTS = np.column_stack((X.ravel(), Y.ravel(), np.zeros(N_X*N_Y)))

lattice_points = EQUILIBRIUM_LATTICE_POINTS.copy()
vx, vy = np.meshgrid(np.zeros(N_X), np.zeros(N_Y))
ion_velocities = np.column_stack((vx.ravel(), vy.ravel()))

electron_2_pos = np.array([- electron_velocity * ELECTRON_2_DELAY, Y_POS_ELECTRON_2, 0.0])
electron_2_vel = np.array([electron_velocity, 0.0, 0.0])
electron_2_coupling = 0.6   # strength of wake interaction

def electron_movement(t):
    return np.array([electron_velocity * t, Y_POS_ELECTRON, 0.0])

def wake_center(electron_1_pos):
    disp = EQUILIBRIUM_LATTICE_POINTS[:, :2] - lattice_points[:, :2]
    mask = lattice_points[:, 0] < electron_1_pos[0]

    if not np.any(mask):
        return electron_1_pos[:2].copy()

    weights = np.linalg.norm(disp[mask], axis=1)
    weighted_positions = EQUILIBRIUM_LATTICE_POINTS[mask, :2] * weights[:, np.newaxis]

    center = weighted_positions.sum(axis=0) / weights.sum()
    return center

def update_electron_2(wake_pos):
    global electron_2_pos, electron_2_vel
    r_vec = wake_pos[:2] - electron_2_pos[:2]
    r2 = r_vec[0]**2 + r_vec[1]**2
    force_vec = electron_2_coupling * r_vec / (regularization + r2**scale_factor)

    electron_2_vel[:2] += force_vec * dt
    electron_2_pos[:2] += electron_2_vel[:2] * dt


def update_ion_lattice(electron_positions):
    global ion_velocities, lattice_points

    for electron_position in electron_positions:
        dxp = lattice_points[:, 0] - electron_position[0]
        dyp = lattice_points[:, 1] - electron_position[1]
        r2 = dxp*dxp + dyp*dyp
        den = regularization + r2**scale_factor

        ion_velocities[:, 0] += -coupling_constant * dxp / den
        ion_velocities[:, 1] += -coupling_constant * dyp / den

    # relaxation (once per step)
    ion_velocities[:, 0] += relaxation_force * (EQUILIBRIUM_LATTICE_POINTS[:,0] - lattice_points[:,0])
    ion_velocities[:, 1] += relaxation_force * (EQUILIBRIUM_LATTICE_POINTS[:,1] - lattice_points[:,1])

    ion_velocities *= 0.98
    lattice_points[:, :2] += ion_velocities * dt



#Plotter setup (pyvista)
plotter = pv.Plotter(off_screen=True, window_size=(1920, 1088))
plotter.set_background((38, 38, 38))
plotter.view_xy(
    bounds=(
        N_X * 0.25, N_X * 0.75,
        N_Y * 0.25, N_Y * 0.75,
        -.5, .5
    )
)


n_setup = n_frames // 10
for i in range(n_setup + 1):
    t = - dt * (n_setup - i)
    update_ion_lattice([electron_movement(t), electron_movement(t - ELECTRON_2_DELAY)])

ion_mesh = pv.PolyData(lattice_points)
glyph_geom = pv.Sphere(radius=0.2)
ion_actor = plotter.add_mesh(
    ion_mesh,
    render_points_as_spheres=True,
    point_size=12,
    color="#84B819",
)


electron_mesh = pv.Sphere(radius=0.1)
electron_actor = plotter.add_mesh(electron_mesh, color="#D10000")

electron_mesh_2 = pv.Sphere(radius=0.1)
electron_actor_2 = plotter.add_mesh(electron_mesh_2, color="#1E90FF")


#Animation loop
with imageio.get_writer(
    "el_ph.mp4",
    fps=fps,
    codec="libx264",
    quality=5
) as writer:

    for i in range(n_frames):
        t = i * dt

        electron_positions = []

        # First electron (driver)
        e1_pos = electron_movement(t)
        electron_actor.SetPosition(*e1_pos)
        electron_positions.append(e1_pos)

        # Compute wake center
        wake_pos = wake_center(e1_pos)

        # Update electron 2 toward wake
        update_electron_2(wake_pos)
        electron_actor_2.SetPosition(*electron_2_pos)
        electron_positions.append(electron_2_pos.copy())


        # Ion dynamics
        update_ion_lattice(electron_positions)

        # Update lattice visualization
        ion_mesh.points[:] = lattice_points
        ion_mesh.Modified()

        # Render
        plotter.render()
        writer.append_data(plotter.screenshot(return_img=True))
