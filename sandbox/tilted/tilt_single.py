import numpy as np
import meep as mp

import matplotlib.pyplot as plt;plt.ion()

############################################
####	User Input ####
############################################

is_vac = [False, True][1]

resolution = 20

rot_angle = np.radians(10)

pol = ['s','p'][0]

dpml = 2
dpad = 3.5
gap = 4
thick = 4

wave = 0.6
df = 0.02  # turn-on bandwidth
n = 1

############################################
############################################

#Derived
fcen = 1/wave
lx_np = dpad*2 + thick
ly_np = dpad*2 + gap
cell_size = mp.Vector3(dpml*2 + lx_np, dpml*2 + ly_np)
src_comp = getattr(mp, {'s': 'Ez', 'p': 'Ey'}[pol])

def pw_amp(k, x0, a):
    def _pw_amp(x):
        return np.exp(1j * k.dot(x + x0)) * a
    return _pw_amp

kdir = mp.Vector3(np.cos(rot_angle), np.sin(rot_angle))  # direction of k (length is irrelevant)
n = 1 # refractive index of material containing the source
k = kdir.unit().scale(2 * np.pi * fcen * n)  # k with correct length

#Source
k_point = mp.Vector3(fcen*n).rotate(mp.Vector3(z=1), rot_angle)
# sim_src = mp.GaussianSource(fcen, fwidth=df, is_integrated=True)
sim_src = mp.ContinuousSource(fcen, fwidth=df, is_integrated=True)
sources = [mp.Source(sim_src, component=src_comp,
            center=mp.Vector3(x=-0.5*lx_np), \
            size=mp.Vector3(y=cell_size.y),
            amp_func=pw_amp(k, mp.Vector3(x=-0.5 * lx_np), kdir.x)),

           mp.Source(sim_src, component=src_comp,
            center=mp.Vector3(y=-0.5*ly_np), \
            size=mp.Vector3(x=cell_size.x),
            amp_func=pw_amp(k, mp.Vector3(y=-0.5 * ly_np), kdir.y)),
        ]

#PML
pml_layers = [mp.PML(thickness=dpml)]#,direction=mp.X)]
# pml_layers = [mp.Absorber(thickness=dpml)]#,direction=mp.X)]

#Geometry
geometry = []
if not is_vac:
    wafer = mp.Block(material=mp.metal, size=mp.Vector3(thick, mp.inf, mp.inf), center=mp.Vector3())
    air = mp.Block(material=mp.air, size=mp.Vector3(thick, gap, mp.inf), center=mp.Vector3())
    geometry += [wafer, air]

#Build simulation
sim = mp.Simulation(force_complex_fields=True,
    resolution=resolution, default_material=mp.Medium(index=n),
    cell_size=cell_size, boundary_layers=pml_layers, sources=sources,
    geometry=geometry, k_point=k_point)

############################################
############################################

#Add DFT fields
vol = mp.Volume(size=mp.Vector3(lx_np, ly_np))
# dft_obj = sim.add_dft_fields([src_comp], [fcen], where=vol)


#Run
# sim.run(until_after_sources=mp.stop_when_fields_decayed(50, src_comp, mp.Vector3(-0.5*cell_size.x+dpml), 1e-6))
sim.run(until=400)
sim.plot2D(fields=src_comp, output_plane=vol)

fld = sim.get_array(src_comp)
fig, axes = plt.subplots(1,2)
axes[0].imshow(np.abs(fld).T, origin='lower')
axes[1].imshow(np.angle(fld).T, origin='lower')

# #Get field
# fld = sim.get_dft_array(dft_obj, src_comp, 0)


breakpoint()
