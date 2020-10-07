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
    # 'sim_geometry':     'edge',
    'wafer_material':   'cSi',
    'wafer_thick':      2.,
    'skin_material':    'Au',
    'skin_thick':       0.4,
    'gap_width':        5,

    # 'taper_angle':      20,
    # 'wall_thick':       0.2,
    # 'scallop_height': 0.8,
    # 'scallop_depth': 0.2,

    ### Numerics ###
    'resolution':       10,
    'padx':             2,
    'pady':             2.,
    'padz':             2.,
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
