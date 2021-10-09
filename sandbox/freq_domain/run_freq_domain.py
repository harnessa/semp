import semp
import meep as mp
import h5py

wave = 0.641
seam_dark = 5
seam_lite = 10
resolution = 20
pml = 4
pad = 4

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

#Load simulation
prop = semp.Propagator(MEEP_params, PROP_params)
sim = prop.msim.build_sim(pol=pol)

#Run
sim.init_sim()
sim.run(until=160)

ez = sim.get_array(component=mp.Ez)

from mpi4py import MPI
rank = MPI.COMM_WORLD.rank

if rank == 0:
    with h5py.File('./edge_time2_ez.h5', 'w') as f:
        f.create_dataset('ez', data=ez)
