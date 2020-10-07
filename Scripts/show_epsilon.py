"""
show_epsilon.py

Author: Anthony Harness
Affiliation: Princeton University
Created on: 10-07-2020
Package: SEMP

Description: Build epsilon geometry for SEMP simulation
License: Refer to $pkg_home_dir/LICENSE
"""

import semp
import numpy as np
import matplotlib.pyplot as plt;plt.ion()

MEEP_params = {
    ### Lab Properties  ###
    'polarization':     ['s','p'][0],
    'wave':             0.641,

    ### Mask Properties ###
    # 'sim_geometry':     'edge',
    'wafer_material':   'cSi',
    'skin_material':    'Au',
    # 'wafer_material':   'metal',
    # 'skin_material':    'metal',
    'wafer_thick':      2.,
    'skin_thick':       0.4,
    'gap_width':        5,

    'taper_angle':      20,
    # 'wall_thick':       0.2,
    # 'scallop_height': 0.8,
    # 'scallop_depth': 0.2,

    'sim_geometry':     'corner',
    'corner_length':    2,
    # 'corner_dy':        0.5,
    # 'corner_dz':        1,

    ### Numerics ###
    'resolution':       10,
    'padx':             1/2,
    'pady':             1/2,
    'padz':             1/2,
    'padx':             1,
    'pady':             1.,
    'padz':             1.,
    'n_periods':        50,
}

PROP_params = {
    'do_save':          True,
    'session_name':     'test',
    'save_ext':         '',
    'is_movie':         True,
}

prop = semp.Propagator(MEEP_params, PROP_params)

eps = np.abs(prop.get_epsilon())
ms = prop.meep_sim


fig, axes = plt.subplots(3,1,figsize=(11,8.5),tight_layout=True)
ff = [0.5,0.5,0.55]
dd = [int(ff[j]*eps.shape[j]) for j in range(3)]
for j in range(3):
    axes[j].imshow(eps.take(dd[j],axis=j), vmax=15, interpolation='none')
axes[0].axhline(dd[1], color='r')
axes[0].axvline(dd[2], color='r')
axes[0].set_xlabel('z')
axes[0].set_ylabel('y')
axes[1].set_xlabel('z')
axes[1].set_ylabel('x')
axes[2].set_xlabel('y')
axes[2].set_ylabel('x')
axes[0].axhline(eps.shape[1]/2, color='k', linestyle=':')
axes[0].axvline(eps.shape[2]/2, color='k', linestyle=':')
axes[1].axhline(eps.shape[0]/2, color='k', linestyle=':')
axes[1].axvline(eps.shape[2]/2, color='k', linestyle=':')
axes[2].axhline(eps.shape[0]/2, color='k', linestyle=':')
axes[2].axvline(eps.shape[1]/2, color='k', linestyle=':')

import pdb;pdb.set_trace()
