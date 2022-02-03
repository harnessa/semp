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

# scallop_list = [
#     [(0.205, 0.324), (0.762, 0.810)],
#     [(0.148, 1.157), (0.762, 0.824)],
#     [(0.0190, 1.900), (0.762, 0.767)]
# ]

MEEP_params = {

    ### Lab Properties  ###
    'waves':            0.641,

    ### Mask Properties ###
    'sim_geometry':     'gap',
    'seam_dark':        1,
    'seam_lite':        1,

    'wafer_material':   'cSi',
    'skin_material':    'Al',
    'wafer_thick':      2.25,
    'skin_thick':       0.3,
    # 'wall_thick':       0.01,

    # 'oxide_material':   'Al2O3',
    # 'oxide_material':   'metal',
    # 'oxide_thick':      0.2,

    # 'scallop_list':     scallop_list,

    ### Numerics ###
    'resolution':       20,
    'pml_all':          2,
    'pad_all':          0,

}

#Load propagator + analyzer
prop = semp.Propagator(MEEP_params, {'verbose':False}, is_analysis=True)
alz = semp.analysis.Analyzer({}, prop=prop, build_geo=True)

#Show epsilon
axes = alz.show_epsilon(with_lines=False)


plt.ioff()
from mpl_toolkits.axes_grid1 import make_axes_locatable
fig, axes = plt.subplots(1, figsize=(8,6))

#Get image extent
extent = [alz.yy[-1], alz.yy[0], alz.xx[-1], alz.xx[0]]
#Adjust for yee lattice
extent = [x + 0.5/alz.prop.msim.resolution for x in extent]

#Show image
vmax = 25
out = axes.imshow(alz.eps, extent=extent, interpolation='none',
    vmin=0, vmax=vmax, cmap=plt.cm.get_cmap('viridis', 3))

#Colorbar
cticks = [0., vmax/2, vmax]
clabs = ['Vacuum', 'Silicon', 'Aluminum']
divider = make_axes_locatable(axes)
cax = divider.append_axes("right", size="3%", pad=0.05)
cbar = plt.colorbar(out,cax=cax, ticks=cticks, orientation='vertical')
cbar.set_ticklabels(clabs)

#Labels
axes.set_xlabel('Y [microns]')
axes.set_ylabel('X [microns]')

axes.axhline(-MEEP_params['pml_all'], color='r', linestyle=':')
axes.axhline(alz.prop.msim.wafer_thick/2, color='k', linestyle=':')
fig.savefig('/home/aharness/Desktop/eps.png')
breakpoint()
