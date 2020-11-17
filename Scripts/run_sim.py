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
    # 'sim_geometry':     'corner',
    # 'wafer_material':   'cSi',
    # 'skin_material':    'Au',
    'wafer_thick':      2.,
    'skin_thick':       0.,
    'gap_width':        3,

    # 'taper_angle':      20,
    # 'wall_thick':       0.2,
    # 'scallop_height': 0.8,
    # 'scallop_depth': 0.2,
    'corner_length':    1,

    ### Numerics ###
    'resolution':       20,
    'padx':             1,
    'pady':             1,
    'padz':             1,
    'n_periods':        40,
    # 'use_absorber':     True,
}

PROP_params = {
    'save_nt':          3,
    'do_save':          True,
    # 'session':          'corner_test',
    'session':          'test',
    'ext':              'full',
    'output_full_dim':  True,
}

prop = semp.Propagator(MEEP_params, PROP_params)
prop.run_sim()
