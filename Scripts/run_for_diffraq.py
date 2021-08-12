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


base_dir = 'Si_2'

seam_dark = 5
seam_lite = 10

wafer_thick = 2
skin_thick = 0.25
skin_metal = 'metal'

resolution = 30
n_periods = 50

waves = [0.641, 0.660, 0.699, 0.725]

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
    PROP_params['session'] = f'{base_dir}/{wv*1e3:.0f}_nm'

    #Run simulation
    prop = semp.Propagator(MEEP_params, PROP_params)
    prop.run_sim()
