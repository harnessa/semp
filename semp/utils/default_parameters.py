"""
default_parameters.py

Author: Anthony Harness
Affiliation: Princeton University
Created on: 10-06-2020
Package: SEMP

Description: Default parameters for the Meep simulation class and the Propagator class.
License: Refer to $pkg_home_dir/LICENSE
"""

import semp
import numpy as np

""" All units in [um] unless specificed others """

#Default base directory
base_dir = r"%s/Results"%semp.pkg_home_dir

#Default parameters for Meep simulation
def_params_MEEP = {

    ### Lab Properties  ###
    'obs_distance':     0.,        # Distance from wafer bottom to near field.
    'polarization':     's',       # Polarization. Options:'s'-skipping, 'p'-plunging
    'wave':             0.5,       # Central wavelength.
    'is_diverging':     False,     # Diverging light source
    'source_distance':  27.5,      # Units: [m], Light source distance for diverging beam
    'source_offset_y':  0,         # Units: [m], Light source center in y
    'source_offset_z':  0,         # Units: [m], Light source center in z

    ### Mask Properties ###
    'sim_geometry':     'gap',     # Geometry of simulation. Options: [gap, edge, corner, vacuum]
    'wafer_material':   'metal',
    'skin_material':    'metal',
    'wafer_epsilon':    None,      # Complex permittivity of wafer material
    'skin_epsilon':     None,      # Complex permittivity of skin material
    'wafer_thick':      1.,
    'skin_thick':       0.,
    'wall_thick':       0.,
    'gap_width':        5.,
    'scallop_height':   0.,        # Height or distance between scallops
    'scallop_depth':    0.,        # Depth of scallops
    'taper_angle':      0,         # Angle to taper edge [deg]
    'is_sommerfeld':    False,     # Is sommerfeld solution model
    'corner_length':    0.,
    'corner_dy':        0.,        # Distance broken corner extrudes in y
    'corner_dz':        0.,        # Distance broken corner extrudes in z

    ### Numerics ###
    'resolution':       30,        # [pixels / um]
    'use_absorber':     False,     # Use absorber instead of PML
    'pmlx':             2.,
    'pmly':             2.,
    'pmlz':             2.,
    'padx':             6.,
    'pady':             6.,
    'padz':             6.,
    'fwid':             0.2,       # Frequency width [multiple of frequency]
    'n_periods':        50,        # Number of optical time periods to run
    'decay_tol':        1e-9,      # Run until decays to this level
    'courant':          0.5,       # Courant number (lower is slower, but more numerically stable)
}

##############################################
##############################################

#Default parameters for Propagator
def_params_PROP = {

    ### Movie ###
    'save_nt':          1,         # Number of saves per optical time period
    'is_movie':         False,     # Run movie?

    ### Saving ###
    'save_dir_base':    semp.results_dir,
    'session_name':     '',
    'save_ext':         '',
    'do_save':          False,     # Save data?
    'save_all_comps':   False,     # Save all field components?

}

##############################################
##############################################

def_params = {'MEEP_params': def_params_MEEP, 'PROP_params': def_params_PROP}
