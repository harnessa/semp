import numpy as np
import semp
import meep as mp
import h5py
import time

wave = 0.641
seam_dark = 5
seam_lite = 10
resolution = 25
pml = 4
pad = 4

run_time = 80

MEEP_params = {
    ### Lab Properties  ###
    'polars':           ['s', 'p'],
    'wave':             wave,

    ### Mask Properties ###
    'sim_geometry':     'edge',
    'gap_width':        5,
    'is_sommerfeld':    True,
    'seam_dark':        seam_dark,
    'seam_lite':        seam_lite,

    ### Numerics ###
    'resolution':       resolution,
    'pml_all':          pml,
    'pad_all':          pad,
    'n_periods':        50,
    'use_absorber':     True,

}

PROP_params = {
    'verbose':          True,
    'base_dir':         f'{semp.tmp_dir}/tests',
    'session':          'sommer',
}

pol = ['s','p'][0]
is_vac = False

#Load simulation
prop = semp.Propagator(MEEP_params, PROP_params)
msim = prop.msim

#K-point
k_point = mp.Vector3()

#PML
pml_layers = msim.get_pml_layers(is_vac)

#Symmetries
symmetries = msim.get_symmetries(is_vac, pol)

#Computational cell
cell_size = msim.get_cell_size(is_vac)

#Source
sources = msim.get_source(pol, is_vac)

#Geometry
geometry = msim.get_geometry(is_vac)

#Build simulation
sim = mp.Simulation(split_chunks_evenly=False, force_complex_fields=True,
    ensure_periodicity=False, resolution=msim.resolution, Courant=msim.courant,
    cell_size=cell_size, boundary_layers=pml_layers, sources=sources,
    geometry=geometry, symmetries=symmetries, k_point=k_point)

nonpml_vol = mp.Volume(mp.Vector3(), size=mp.Vector3(msim.geo.lx-2*msim.geo.pmlx, msim.geo.ly-2*msim.geo.pmly))

##########################################
tik = time.perf_counter()

#Normal run
sim.init_sim()
sim.run(until=run_time)

ez = sim.get_array(vol=nonpml_vol, component=mp.Ez)

from mpi4py import MPI
rank = MPI.COMM_WORLD.rank

if rank == 0:
    print(f'Time 1: {time.perf_counter()-tik:.3f}')
    
    with h5py.File('./normal_ez.h5', 'w') as f:
        f.create_dataset('ez', data=ez)

###########################################
tik = time.perf_counter()

#Gaussian run
sim.reset_meep()

#Center of source
src_pt = mp.Vector3(x=msim.geo.source_x)

#Size of source
src_sze_y = [msim.geo.ly, 0.][is_vac]
src_sze_z = [msim.geo.lz, 0.][is_vac]

#Get source functions
sim_src = mp.GaussianSource(msim.fcen, fwidth=0.1, is_integrated=True)

#Get source component
src_comp = getattr(mp, {'s': 'Ez', 'p': 'Hz'}[pol])

#Build source   #TODO: add gaussian beam source option
sources = [mp.Source(sim_src, component=src_comp, center=src_pt, \
    size=mp.Vector3(y=src_sze_y, z=src_sze_z), amp_func=None)]

sim = mp.Simulation(split_chunks_evenly=False, force_complex_fields=True,
    ensure_periodicity=False, resolution=msim.resolution, Courant=msim.courant,
    cell_size=cell_size, boundary_layers=pml_layers, sources=sources,
    geometry=geometry, symmetries=symmetries, k_point=k_point)

dft_obj = sim.add_dft_fields([mp.Ez], msim.fcen, 0, 1, where=nonpml_vol)

sim.init_sim()
sim.run(until=run_time)

ez2 = sim.get_dft_array(dft_obj, mp.Ez, 0)

from mpi4py import MPI
rank = MPI.COMM_WORLD.rank

if rank == 0:
    print(f'Time 2: {time.perf_counter()-tik:.3f}')

    with h5py.File('./dft_ez.h5', 'w') as f:
        f.create_dataset('ez', data=ez2)

###########################################
