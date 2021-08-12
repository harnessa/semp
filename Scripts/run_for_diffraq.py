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


base_dir = 'comp_M12P3'

resolution = 50

wafer_thick = 2
skin_metal = 'Al'
skin_thick = 0.25

scallop_height = 0.8
scallop_depth  = 0.225
taper_angle = 5

seam_dark = 10
seam_lite = 20
n_periods = 150

waves = [0.641, 0.660, 0.699, 0.725][:1]

######################################

#Meep parameters
MEEP_params = {

    ### Lab Properties  ###
    'polars':           ['s', 'p'],

    ### Mask Properties ###
    'sim_geometry':     'edge',
    'seam_dark':        seam_dark,
    'seam_lite':        seam_lite,

    'wafer_material':   'cSi',
    'skin_material':    skin_metal,
    'wafer_thick':      wafer_thick,
    'skin_thick':       skin_thick,

    'scallop_height':   scallop_height,
    'scallop_depth':    scallop_depth,
    'taper_angle':      taper_angle,

    ### Numerics ###
    'resolution':       resolution,
    'pml_all':          4,
    'pad_all':          4,
    'n_periods':        n_periods,

}

#Main parameters
PROP_params = {
}

#Loop through wavelengths and run
for wv in waves:

    #Add wavelength to parameters
    MEEP_params['wave'] = wv
    PROP_params['session'] = f'{base_dir}/{wv*1e3:.0f}nm'

    #Run simulation
    prop = semp.Propagator(MEEP_params, PROP_params)
    prop.run_sim()