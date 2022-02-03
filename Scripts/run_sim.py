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
    'polars':           ['s', 'p'][0],
    'waves':            0.641,

    ### Mask Properties ###
    'sim_geometry':     'gap',
    'seam_dark':        4,
    'seam_lite':        4,
    'is_sommerfeld':    False,

    # 'wafer_material':   'cSi',
    'wafer_epsilon':    (14.960 + 1j*0.14014),
    # 'wafer_epsilon':    (14.960 + 1j*0),
    # 'skin_epsilon':    (-44.4 + 1j*16),
    'skin_material':    'Al',
    'wafer_thick':      2,
    'skin_thick':       0.25,

    ### Numerics ###
    'resolution':       30,
    'pml_all':          2,
    'pad_all':          2,
    'use_absorber': False,
}

#Main parameters
PROP_params = {
    # 'session': 'sommer'
    'session':  'skin_epsi_2',
    # 'session': 'gauss_beam/edge_gg_off',
    # 'is_movie':True,
}

#Run simulation
prop = semp.Propagator(MEEP_params, PROP_params)
prop.run_sim()
