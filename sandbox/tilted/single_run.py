import numpy as np
import meep as mp

import matplotlib.pyplot as plt;plt.ion()

############################################
####	User Input ####
############################################

is_vac = [False, True][1]

resolution = 20

rot_angle = np.radians(0)

pol = ['s','p'][0]

dpml = 1
dpad = 1
dgap = 3
thick = 1

wave = 0.6
df = 0.02  # turn-on bandwidth
n = 1

############################################
############################################

#Derived
fcen = 1/wave
cell_size = mp.Vector3(dpml*2 + dpad*2 + thick, dpml*2 + dpad*2 + dgap)
src_comp = getattr(mp, {'s': 'Ez', 'p': 'Ey'}[pol])

#Source
# k_point = mp.Vector3(fcen*n).rotate(mp.Vector3(z=1), rot_angle)
k_point = mp.Vector3()
sim_src = mp.GaussianSource(fcen, fwidth=df, is_integrated=True)
# sim_src = mp.ContinuousSource(fcen, fwidth=df, is_integrated=True)
sources = [mp.Source(sim_src, component=src_comp, center=mp.Vector3(x=-0.5*cell_size.x+dpml), \
    size=mp.Vector3(y=cell_size.y), amp_func=None)]

# sources = [mp.EigenModeSource(src=sim_src,
#                               center=mp.Vector3(x=-0.5*cell_size.x+dpml),
#                               size=mp.Vector3(y=cell_size.y),
#                               direction=mp.AUTOMATIC if rot_angle == 0 else mp.NO_DIRECTION,
#                               eig_kpoint=k_point,
#                               eig_band=1,
#                               eig_parity=mp.EVEN_Y+mp.ODD_Z if rot_angle == 0 else mp.ODD_Z,
#                               eig_match_freq=True)]


#PML
pml_layers = [mp.PML(thickness=dpml)]#,direction=mp.X)]

#Geometry
geometry = []
if not is_vac:
    wafer = mp.Block(material=mp.metal, size=mp.Vector3(thick, mp.inf, mp.inf), center=mp.Vector3())
    air = mp.Block(material=mp.air, size=mp.Vector3(thick, dgap, mp.inf), center=mp.Vector3())
    geometry += [wafer, air]

#Build simulation
sim = mp.Simulation(force_complex_fields=True,
    ensure_periodicity=False, resolution=resolution,
    cell_size=cell_size, boundary_layers=pml_layers, sources=sources,
    geometry=geometry, symmetries=[], k_point=k_point)

############################################
############################################

#Add DFT fields
vol = mp.Volume(size=mp.Vector3(dpad*2 + thick, dpad*2 + dgap))
dft_obj = sim.add_dft_fields([src_comp], [fcen], where=vol)

#Set output directory + prefix
import os
sim.use_output_directory(os.getcwd())
sim.filename_prefix = ''

#Get filename
fld_name = {'s':'ez', 'p':'hz'}[pol]
filename = f'movie_{fld_name}'
#Function to save field (appended) to h5 file
if pol == 's':
    fld_func = mp.output_efield_z
else:
    fld_func = mp.output_hfield_z
h5_func = mp.to_appended(filename, mp.at_every(1, fld_func))

#Run
sim.run(h5_func, until_after_sources=mp.stop_when_fields_decayed(50, src_comp, mp.Vector3(-0.5*cell_size.x+dpml), 1e-6))
# sim.run(h5_func, until=100)

#Get field
fld = sim.get_dft_array(dft_obj, src_comp, 0)
x,y,z,w = sim.get_array_metadata(dft_cell=dft_obj)

sim.plot2D(fields=src_comp, output_plane=vol)
fig, axes = plt.subplots(1,2)
axes[0].imshow(np.abs(fld).T)
axes[1].imshow(np.angle(fld).T)

dd = np.angle(fld)[:,50] % (2*np.pi)
ee = np.angle(np.exp(1j*2*np.pi*fcen*x))
ee += dd[5] - ee[5]
ee %= 2*np.pi
plt.figure()
plt.plot(x, dd)
plt.plot(x, ee, '--')

import h5py
with h5py.File(filename+'.h5', 'r') as f:
    et = f['ez.r'][()] + 1j*f['ez.i'][()]

tf = np.angle(et[et.shape[0]//2, et.shape[1]//2])
tf %= 2*np.pi
tt = np.arange(len(tf)) * 0.5*resolution
tf2 = np.angle(np.exp(-1j*2*np.pi*fcen*tt))
tf2 += tf[15] - tf2[15]
tf2 %= 2*np.pi

plt.figure()
plt.plot(tt, tf)
plt.plot(tt, tf2, '--')

breakpoint()
