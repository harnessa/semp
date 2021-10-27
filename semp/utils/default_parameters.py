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
base_dir = rf"{semp.pkg_home_dir}/Results"

############################################
####	MEEP Parameters ####
############################################

#Default parameters for Meep simulation
def_params_MEEP = {

    ### Lab Properties  ###
    'polars':           ['s','p'],  # Polarization. Options:'s'-skipping, 'p'-plunging
    'waves':            np.array([0.5]),        # Array of wavelengths.
    'is_diverging':     False,      # Diverging light source
    'source_distance':  27.5,       # Units: [m], Light source distance for diverging beam
    'source_offset_y':  0,          # Units: [m], Light source center in y
    'source_offset_z':  0,          # Units: [m], Light source center in z

    ### Mask Properties ###
    'sim_geometry':     'edge',     # Geometry of simulation. Options: [gap, edge, corner, vacuum]
    'wafer_material':   'metal',    # Wafer material
    'skin_material':    'metal',    # Skin and sidewall material
    'oxide_material':   'metal',    # Oxidation material
    'wafer_epsilon':    None,       # Complex permittivity of wafer material
    'skin_epsilon':     None,       # Complex permittivity of skin material
    'oxide_epsilon':    None,       # Complex permittivity of oxide material
    'is_sommerfeld':    False,      # Is sommerfeld solution model
    'wafer_thick':      1.,         # Wafer thickness [um]
    'skin_thick':       0.,         # Skin thickness [um]
    'oxide_thick':      0.,         # Oxidation thickness [um]
    'wall_thick':       0.,         # Sidewall thickness [um]
    'gap_width':        5.,         # Width of gap (also w/ edge) [um]
    'seam_dark':        5.,         # Width of seam on dark side (wafer) [um]
    'seam_lite':        5.,         # Width of seam on light side (vacuum) [um]
    'scallop_height':   0.,         # Height or distance between scallops [um]
    'scallop_depth':    0.,         # Depth of scallops (>0 turns on scallops) [um]
    'scallop_start':    0.,         # Vertical distance below skin where scallops start [um]
    'taper_angle':      0,          # Angle to taper edge [deg]
    'shave_angle':      0,          # Angle to shave off scallops [deg]
    'footing_size':     None,       # Height and width of ellipsoidal footing [um]
    'corner_length':    0.,         # Length of gap extending in z  [um]
    'corner_dy':        0.,         # Distance broken corner extrudes in y [um]
    'corner_dz':        0.,         # Distance broken corner extrudes in z [um]
    'scallop_list':     [],         # List of scallops. Ellipsoid center, ellipsoid size [um]
    'scallop_ball':     0,          # Radius of ball on scallop [um]

    ### Numerics ###
    'resolution':       30,         # [pixels / um]
    'use_absorber':     True,       # Use absorber instead of PML
    'pmlx':             2.,         # Size of PML in X [um]
    'pmly':             2.,         # Size of PML in Y [um]
    'pmlz':             2.,         # Size of PML in Z [um]
    'padx':             6.,         # Vertical padding between PML and wafer [um]
    'pady':             6.,         # Horizontal padding between PML and gap [um]
    'padz':             6.,         # 3D padding between PML and start of gap [um]
    'pml_all':          None,       # If not None, replaces all PML components with value
    'pad_all':          None,       # If not None, replaces all pad components with value
    'n_periods':        50,         # Number of optical time periods to run
    'decay_dt':         50,         # Time after decay to
    'decay_by':         1e-3,       # Decay amount
    'courant':          0.5,        # Courant number (lower is slower, but more numerically stable)
}

##############################################
##############################################

############################################
####	PROP Parameters ####
############################################

#Default parameters for Propagator
def_params_PROP = {

    ### Movie ###
    'save_nt':          1,          # Number of saves per optical time period
    'is_movie':         False,      # Run movie?

    ### Saving ###
    'base_dir':         semp.results_dir,       # Directory base
    'session':          '',         # Session: save under 'base_dir/session'
    'verbose':          True,       # Print statements?
    'save_all':         True,
}

##############################################
##############################################

############################################
####	ANLZ Parameters ####
############################################

#Default parameters for Analyzer
def_params_ANLZ = {
    ### Saving ###
    'do_save':          False,      # Save analysis results?
    'save_dir':         None,       # If None, save to load dir
    'save_ext':         '',         # Save analysis extension
    ### Loading ###
    'base_dir':         semp.results_dir,       # Directory base
    'session':          '',         # Session: load from 'base_dir/session'
    'time_ext':         None,
    ### Analyzing ###
    'obs_distance':     0.,         # Distance from wafer bottom to near field.
}

##############################################
##############################################

def_params = {'MEEP_params': def_params_MEEP, 'PROP_params': def_params_PROP}
