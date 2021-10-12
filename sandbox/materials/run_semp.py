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
    'resolution':       35,
    'pml_all':          4,
    'pad_all':          4,

}

# #Meep material
# prop = semp.Propagator(MEEP_params, {'session':'materials/meep'})
# prop.run_sim()
#
# nn = 3.859
# kk = 0.015
# eps = (nn**2 - kk**2) + 1j*(2*nn*kk)
#
# #Epsilon
# MEEP_params['wafer_epsilon'] = eps
# prop = semp.Propagator(MEEP_params, {'session':'materials/epsilon'})
# prop.run_sim()

#skin Epsilon
# snn =
# skk = 0.015
# seps = (snn**2 - skk**2) + 1j*(2*snn*skk)
seps = -52.55990102266269+20.335479148581857j

MEEP_params['skin_epsilon'] = seps
prop = semp.Propagator(MEEP_params, {'session':'materials/skin_meep_epsilon'})
prop.run_sim()


# #Meep material
# MEEP_params['wafer_material'] = 'mSi'
# prop = semp.Propagator(MEEP_params, {'session':'materials/msi_300Al'})
# prop.run_sim()
