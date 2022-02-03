import numpy as np
import meep as mp
import meep.materials as mat_lib
import h5py

import matplotlib.pyplot as plt;plt.ion()

############################################
####	User Input ####
############################################

is_vac = [False, True][0]

resolution = 30

pol = ['s','p'][0]

dpml = 2
dpad = 3
gap = 1
thick = 5

wave = 0.6
df = 0.1  # turn-on bandwidth

waf_mat = mat_lib.cSi

############################################
############################################

#Derived
fcen = 1/wave
lx_np = dpad*2 + thick
ly_np = dpad*2 + gap
cell_size = mp.Vector3(dpml*2 + lx_np, dpml*2 + ly_np)
src_comp = getattr(mp, {'s': 'Ez', 'p': 'Ey'}[pol])

#Source
sim_src = mp.GaussianSource(fcen, fwidth=df, is_integrated=True)
# sim_src = mp.ContinuousSource(fcen, fwidth=df, is_integrated=True)
src_pt = mp.Vector3(x=-0.5*lx_np)
sources = [mp.Source(sim_src, component=src_comp,
            center=src_pt, \
            size=mp.Vector3(y=cell_size.y))]

#PML
pml_layers = [mp.PML(thickness=dpml)]

#Geometry
geometry = []
if not is_vac:
    met_thick = 0.#25
    wafer2 = mp.Block(material=mp.metal, size=mp.Vector3(met_thick, mp.inf, mp.inf),
        center=mp.Vector3(0.5*cell_size.x - dpml - met_thick/2))
    wafer = mp.Block(material=waf_mat, size=mp.Vector3(thick, mp.inf, mp.inf),
        center=mp.Vector3(0.5*cell_size.x - dpml - met_thick - thick/2))
    geometry += [wafer]#,wafer2]

#Build simulation
sim = mp.Simulation(force_complex_fields=True,
    resolution=resolution, ensure_periodicity=False,
    cell_size=cell_size, boundary_layers=pml_layers, sources=sources,
    geometry=geometry, k_point=mp.Vector3())

# sim.init_sim()
# eps = sim.get_epsilon(fcen)
# plt.imshow(abs(eps),vmax=20)
# breakpoint()
############################################
############################################

#Add DFT fields
vol = mp.Volume(size=mp.Vector3(lx_np, ly_np))
all_flds = [mp.Ex, mp.Ey, mp.Ez, mp.Hx, mp.Hy, mp.Hz]
dft_obj = sim.add_dft_fields(all_flds, [fcen], where=vol)

flux_box = sim.add_flux(fcen, 0, 1,
    mp.FluxRegion(center=mp.Vector3(0.5*cell_size.x - dpml - 0.1),
    size=mp.Vector3(0,ly_np),weight=+1))

#Run
sim.run(until_after_sources=mp.stop_when_fields_decayed(30, src_comp, \
    mp.Vector3(-0.5*cell_size.x+dpml), 1e-5))
# sim.run(until=200)
# sim.plot2D(mp.Ez)

Dz = sim.get_dft_array(dft_obj,mp.Dz,0)
Ez = sim.get_dft_array(dft_obj,mp.Ez,0)
absorbed_power_density = 2*np.pi*fcen * np.imag(np.conj(Ez)*Dz)
absorbed_flux = mp.get_fluxes(flux_box)[0]


# fld = sim.get_array(src_comp, vol=vol)
fld = sim.get_dft_array(dft_obj, src_comp, 0)
x,y,z,w = sim.get_array_metadata(dft_cell=dft_obj)

plt.imshow(abs(fld))
# plt.figure()
# plt.imshow(absorbed_power_density)

breakpoint()
