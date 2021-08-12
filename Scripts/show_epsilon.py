"""
show_epsilon.py

Author: Anthony Harness
Affiliation: Princeton University
Created on: 08-12-2021
Package: SEMP

Description: Build simulation and show dielectric
License: Refer to $pkg_home_dir/LICENSE
"""

import numpy as np
import semp
import matplotlib.pyplot as plt;plt.ion()

MEEP_params = {

    ### Lab Properties  ###
    'wave':             0.641,

    ### Mask Properties ###
    'sim_geometry':     'edge',
    'seam_dark':        5,
    'seam_lite':        10,

    'wafer_material':   'cSi',
    'skin_material':    'metal',
    'wafer_thick':      2,
    'skin_thick':       0.25,

    ### Numerics ###
    'resolution':       30,
    'pml_all':          4,
    'pad_all':          4,

}

#Load propagator + analyzer
prop = semp.Propagator(MEEP_params, {'verbose':False}, is_analysis=True)
alz = semp.analysis.Analyzer({}, prop=prop, build_geo=True)

#Show epsilon
alz.show_epsilon()

breakpoint()
