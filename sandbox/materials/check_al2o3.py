import numpy as np
import semp
import h5py

#Meep parameters
MEEP_params = {
    ### Lab Properties  ###
    'polars':           ['s', 'p'],
    'wave':             0.641,

    ### Mask Properties ###
    'sim_geometry':     'edge',
    'seam_dark':        5,
    'seam_lite':        10,

    'wafer_material':   'cSi',
    'skin_material':    'Al',
    'wafer_thick':      2,
    'skin_thick':       0.25,

    ### Numerics ###
    'resolution':       50,
    'pml_all':          4,
    'pad_all':          4,

}

#Aluminum
MEEP_params['skin_material'] = 'Al'
ext1 = 'al'
prop = semp.Propagator(MEEP_params, {'session':f'materials/{ext1}'})
# prop.run_sim()

#Aluminum oxide
MEEP_params['skin_material'] = 'Al2O3'
ext2 = 'al2o3'
prop = semp.Propagator(MEEP_params, {'session':f'materials/{ext2}'})
# prop.run_sim()

#Load and save
if semp.mpi_rank == 0:

    #Load analyzers
    alz1 = semp.analysis.Analyzer({'session':f'materials/{ext1}'})
    alz2 = semp.analysis.Analyzer({'session':f'materials/{ext2}'})

    #Get data
    is_bbek = True

    with h5py.File(f'./saves/{ext1}.h5', 'w') as f:
        f.create_dataset('yy', data=alz1.yy)
        f.create_dataset('xind', data=alz1.get_xind(obs_x=alz1.prop.msim.wafer_thick/2))
        f.create_dataset('ez', data=alz1.get_data('ez', is_bbek=is_bbek))
        f.create_dataset('hz', data=alz1.get_data('hz', is_bbek=is_bbek))

    with h5py.File(f'./saves/{ext2}.h5', 'w') as f:
        f.create_dataset('yy', data=alz2.yy)
        f.create_dataset('xind', data=alz2.get_xind(obs_x=alz2.prop.msim.wafer_thick/2))
        f.create_dataset('ez', data=alz2.get_data('ez', is_bbek=is_bbek))
        f.create_dataset('hz', data=alz2.get_data('hz', is_bbek=is_bbek))
