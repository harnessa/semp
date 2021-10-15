import numpy as np
import semp
import meep as mp
import h5py
import time
from mpi4py import MPI
rank = MPI.COMM_WORLD.rank

seam_dark = 5
seam_lite = 10
resolution = 30
pml = 12
pad = 12

waves = np.array([0.641, 0.660, 0.699, 0.725])
fcens = 1/waves

run_time = 80

nn = 3.859
kk = 0.015
eps = (nn**2 - kk**2) + 1j*(2*nn*kk)

MEEP_params = {
    ### Mask Properties ###
    'sim_geometry':     'edge',
    'seam_dark':        seam_dark,
    'seam_lite':        seam_lite,

    # 'is_sommerfeld':    True,
    'wafer_material':   'cSi',
    # 'wafer_epsilon':    eps,
    'skin_material':    'Al',
    'wafer_thick':      2,
    'skin_thick':       0.25,


    ### Numerics ###
    'resolution':       resolution,
    'pml_all':          pml,
    'pad_all':          pad,
    'use_absorber':     True,

}

PROP_params = {
    'verbose':          True,
}

pol = ['s','p'][0]
is_vac = False

#Load simulation
prop = semp.Propagator(MEEP_params, PROP_params)
msim = prop.msim

nonpml_vol = mp.Volume(mp.Vector3(), size=mp.Vector3(msim.geo.lx-2*msim.geo.pmlx, msim.geo.ly-2*msim.geo.pmly))

###########################################

def get_sim(sim_src):

    #K-point
    k_point = mp.Vector3()

    #PML
    pml_layers = msim.get_pml_layers(is_vac)

    #Symmetries
    symmetries = msim.get_symmetries(is_vac, pol)

    #Computational cell
    cell_size = msim.get_cell_size(is_vac)

    #Center of source
    src_pt = mp.Vector3(x=msim.geo.source_x)

    #Size of source
    src_sze_y = [msim.geo.ly, 0.][is_vac]
    src_sze_z = [msim.geo.lz, 0.][is_vac]

    #Get source component
    src_comp = getattr(mp, {'s': 'Ez', 'p': 'Hz'}[pol])

    #Build source   #TODO: add gaussian beam source option
    sources = [mp.Source(sim_src, component=src_comp, center=src_pt, \
        size=mp.Vector3(y=src_sze_y, z=src_sze_z), amp_func=None)]

    #Geometry
    geometry = msim.get_geometry(is_vac)

    #Build simulation
    sim = mp.Simulation(split_chunks_evenly=False, force_complex_fields=True,
        ensure_periodicity=False, resolution=msim.resolution, Courant=msim.courant,
        cell_size=cell_size, boundary_layers=pml_layers, sources=sources,
        geometry=geometry, symmetries=symmetries, k_point=k_point)

    return sim

###########################################
###########################################
#Run narrow

# tik = time.perf_counter()
#
# all_ez = []
# for iw in range(len(waves))[:1]:
#
#     sim_src = mp.ContinuousSource(fcens[iw], is_integrated=True)
#
#     sim = get_sim(sim_src)
#     sim.run(until=run_time)
#
#     ez = sim.get_array(vol=nonpml_vol, component=mp.Ez)
#
#     all_ez.append(ez)
#
# MPI.COMM_WORLD.Barrier()
#
# if rank == 0:
#     print(f'Time 1: {time.perf_counter()-tik:.3f}')
#
#     with h5py.File('./saves/si_narrow_ez_eps.h5', 'w') as f:
#         f.create_dataset('waves', data=waves)
#         f.create_dataset('ez', data=all_ez)

###########################################
###########################################
#Run broad

tik = time.perf_counter()

sim_src = mp.GaussianSource(fcens.mean(), fwidth=fcens.ptp()*2, is_integrated=True)
sim = get_sim(sim_src)
dft_obj = sim.add_dft_fields([mp.Ez], fcens, where=nonpml_vol)
# sim.run(until=run_time)
pt = mp.Vector3(x=1)
sim.run(until_after_sources=mp.stop_when_fields_decayed(dt=50, pt=pt, c=mp.Ez, decay_by=1e-3))
all_ez = [sim.get_dft_array(dft_obj, mp.Ez, i) for i in range(len(waves))]

MPI.COMM_WORLD.Barrier()

if rank == 0:
    print(f'Time 2: {time.perf_counter()-tik:.3f}')

    with h5py.File('./saves/si_broad_ez_pml12.h5', 'w') as f:
        f.create_dataset('waves', data=waves)
        f.create_dataset('ez', data=all_ez)
