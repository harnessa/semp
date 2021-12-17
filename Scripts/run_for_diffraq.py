"""
run_for_diffraq.py

Author: Anthony Harness
Affiliation: Princeton University
Created on: 08-12-2021
Package: SEMP

Description: Run SEMP simulations at all wavelengths for to be used for DIFFRAQ
License: Refer to $pkg_home_dir/LICENSE
"""

import numpy as np
import semp


session = 'test'

resolution = 20

wafer_thick = 2.25
skin_metal = 'Al'
skin_thick = 0.3

scallop_height = 0.6
scallop_depth  = 0.1
taper_angle = 0

seam_dark = 15
seam_lite = 35

waves = [0.641, 0.660, 0.699, 0.725]

######################################

#Meep parameters
MEEP_params = {

    ### Lab Properties  ###
    'polars':           ['s', 'p'],
    'waves':            waves,

    ### Mask Properties ###
    'sim_geometry':     'edge',
    'seam_dark':        seam_dark,
    'seam_lite':        seam_lite,

    'wafer_material':   'metal',
    'skin_material':    skin_metal,
    'wafer_thick':      wafer_thick,
    'skin_thick':       skin_thick,

    # 'scallop_height':   scallop_height,
    # 'scallop_depth':    scallop_depth,
    # 'taper_angle':      taper_angle,

    ### Numerics ###
    'resolution':       resolution,
    'pml_all':          8,
    'pad_all':          8,

}

#Main parameters
PROP_params = {
    'session':          session,
    'save_all':         False,
}

#Run simulation
prop = semp.Propagator(MEEP_params, PROP_params)
prop.run_sim()
