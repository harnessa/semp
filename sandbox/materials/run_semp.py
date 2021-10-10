import numpy as np
import semp

#Meep parameters
MEEP_params = {
    ### Lab Properties  ###
    'polars':           ['s', 'p'],
    'waves':             0.641,

    ### Mask Properties ###
    'sim_geometry':     'edge',
    'seam_dark':        5,
    'seam_lite':        10,

    'wafer_material':   'cSi',
    'skin_material':    'Al',
    'wafer_thick':      2,
    'skin_thick':       0.25,

    ### Numerics ###
    'resolution':       30,
    'pml_all':          4,
    'pad_all':          4,

}

#Meep material
prop = semp.Propagator(MEEP_params, {'session':'materials/meep'})
prop.run_sim()

nn = 3.859
kk = 0.015
eps = (nn**2 - kk**2) + 1j*(2*nn*kk)

#Epsilon
MEEP_params['wafer_epsilon'] = eps
prop = semp.Propagator(MEEP_params, {'session':'materials/epsilon'})
prop.run_sim()
