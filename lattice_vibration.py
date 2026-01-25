import numpy as np
import pyvista as pv
import imageio

# Lattice setup
N = 16
k = [np.pi / 4, np.pi / 8]

X, Y = np.meshgrid(np.arange(N), np.arange(N))
initial_lattice_points = np.column_stack((X.ravel(), Y.ravel(), np.zeros(N*N)))

#Video parameters
video_length = 10.0   # seconds
fps = 20              # frames per second
n_frames = int(video_length * fps)
dt = 1 / fps           # time step per frame
num_rotations = 4     # spins precess 20 full turns over the video
omega = 2 * np.pi * num_rotations / video_length

displacement = 0.2

def lattice_positions(t):
    phase = k[0] * X + k[1] * Y - omega * t
    dx = displacement * np.sin(phase)
    return np.column_stack(( (X + dx).ravel(), (Y + dx).ravel(), np.zeros(N*N) ))

#Plotter setup (pyvista)
plotter = pv.Plotter(off_screen=True, window_size=(1920, 1088))
plotter.set_background((38, 38, 38))
plotter.view_xy(bounds=(0, N, 0, N, -1, 1))


mesh = pv.PolyData(lattice_positions(0))
ions = mesh.glyph(
    geom=pv.Sphere(radius=0.2),
    scale=False,
    orient=False
)
actor = plotter.add_mesh(ions, color="#84B819", smooth_shading=True)


#Animation loop
frames = []

for i in range(n_frames):
    t = i * dt
    
    mesh = pv.PolyData(lattice_positions(t))
    ions = mesh.glyph(
        geom=pv.Sphere(radius=0.2),
        scale=False,
        orient=False
    )
    
    actor.mapper.SetInputData(ions)
    plotter.render()
    frames.append(plotter.screenshot(return_img=True))

plotter.close()

#Save Video
imageio.mimsave(
    "lattice_vibration.mp4",
    frames,
    fps=fps,
    codec="libx264",
    quality=5
)
