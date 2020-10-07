"""
run_sim.py

Author: Anthony Harness
Affiliation: Princeton University
Created on: 10-06-2020
Package: SEMP

Description: Run SEMP simulation
License: Refer to $pkg_home_dir/LICENSE
"""

import semp
import numpy as np

MEEP_params = {
    ### Lab Properties  ###
    'polarization':     ['s','p'][0],
    'wave':             0.641,
    ### Mask Properties ###
    'wafer_material':   'cSi',
    'wafer_thick':      3.,
    'skin_material':    'Au',
    'skin_thick':       0.4,
    ### Numerics ###
    'resolution':       30,
    'pmlx':             2,
    'padx':             6,
    'n_periods':        50,
}

PROP_params = {
    'do_save':          True,
    'session_name':     'test',
    'save_ext':         '',
    'is_movie':         True,
}

prop = semp.Propagator(MEEP_params, PROP_params)
prop.run_sim()
