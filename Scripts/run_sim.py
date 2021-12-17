"""
run_sim.py

Author: Anthony Harness
Affiliation: Princeton University
Created on: 08-11-2021
Package: SEMP

Description: Run SEMP simulation
License: Refer to $pkg_home_dir/LICENSE
"""

import numpy as np
import semp

#Meep parameters
MEEP_params = {
    ### Lab Properties  ###
    'polars':           ['s', 'p'],
    'waves':            0.641,

    ### Mask Properties ###
    'sim_geometry':     'edge',
    'seam_dark':        5,
    'seam_lite':        5,
    'is_sommerfeld':    True,

    'wafer_material':   'metal',
    'skin_material':    'metal',
    'wafer_thick':      2,
    # 'skin_thick':       0.25,

    ### Numerics ###
    'resolution':       26,
    'pml_all':          2,
    'pad_all':          2,
    'use_absorber': False,
}

#Main parameters
PROP_params = {
    # 'session': 'sommer'
    'session':  'dum',
    # 'session': 'gauss_beam/edge_gg_off',
    # 'is_movie':True,
}

#Run simulation
prop = semp.Propagator(MEEP_params, PROP_params)
prop.run_sim()
