import meep as mp
import numpy as np
from matplotlib import pyplot as plt; plt.ion()

resolution = 10 # pixels/μm

sx = 14
sy = 10
sz = 0
cell_size = mp.Vector3(sx,sy,sz)

pml_layers = [mp.PML(thickness=2)]

# rotation angle (in degrees) of planewave, counter clockwise (CCW) around z-axis
rot_angle = np.radians(0)

fsrc = 1.0 # frequency of planewave (wavelength = 1/fsrc)

n = 1.5  # refractive index of homogeneous material
default_material = mp.Medium(index=n)

k_point = mp.Vector3(fsrc*n).rotate(mp.Vector3(z=1), rot_angle)

# beam waist
w = 0.6

def amp_func(pt):
    return 1/ np.sqrt(np.pi*w**2) * np.exp(-pt.y ** 2 / w ** 2)

sources = [mp.EigenModeSource(src=mp.GaussianSource(fsrc,fwidth=0.1*fsrc),
                              center=mp.Vector3(),
                              size=mp.Vector3(y=sy),
                              direction=mp.NO_DIRECTION,
                              eig_kpoint=k_point,
                              eig_band=1,
                              eig_parity=mp.EVEN_Y+mp.ODD_Z,
                              amp_func=amp_func,
                              eig_match_freq=True)]

sim = mp.Simulation(cell_size=cell_size,
                    resolution=resolution,
                    boundary_layers=pml_layers,
                    sources=sources,
                    default_material=default_material)

sim.run(until=60)

# plot the intenisty
f = lambda x: np.abs(x) ** 2
field_parameters = {'post_process':f}
sim.plot2D(fields=mp.Ez,field_parameters=field_parameters)
plt.show()

breakpoint()
