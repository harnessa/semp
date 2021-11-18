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

scallop_list = [
    [(0.205, 0.324), (0.762, 0.810)],
    [(0.148, 1.157), (0.762, 0.824)],
    [(0.0190, 1.900), (0.762, 0.767)]
]

MEEP_params = {

    ### Lab Properties  ###
    'waves':            0.641,

    ### Mask Properties ###
    'sim_geometry':     'edge',
    'seam_dark':        15,
    'seam_lite':        35,

    'wafer_material':   'cSi',
    'skin_material':    'Al',
    'wafer_thick':      2.25,
    'skin_thick':       0.3,
    'wall_thick':       0.01,

    'scallop_list':     scallop_list,

    # 'scallop_ball':     0.1,

    ### Numerics ###
    'resolution':       20,
    'pml_all':          8,
    'pad_all':          8,

}

#Load propagator + analyzer
prop = semp.Propagator(MEEP_params, {'verbose':False}, is_analysis=True)
alz = semp.analysis.Analyzer({}, prop=prop, build_geo=True)

#Show epsilon
axes = alz.show_epsilon(with_lines=False)

breakpoint()
