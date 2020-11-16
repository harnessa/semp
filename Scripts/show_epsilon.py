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
    # 'wafer_material':   'cSi',
    # 'skin_material':    'Al',
    'wafer_material':   'metal',
    'skin_material':    'metal',
    'wafer_thick':      2.,
    'skin_thick':       0.4,
    'gap_width':        5,

    # 'taper_angle':      20,
    # 'wall_thick':       0.4,
    # 'scallop_height':   0.8,
    # 'scallop_depth':    0.2,

    'sim_geometry':     'corner',
    'corner_length':    4,
    'corner_dy':        1.5,
    'corner_dz':        1.5,

    ### Numerics ###
    'resolution':       10,
    'pmlx':             1/2,
    'pmly':             1/2,
    'pmlz':             1/2,
    'padx':             1/2,
    'pady':             1/2,
    'padz':             1/2,
}

PROP_params = {
    'do_save':          False,
}

prop = semp.Propagator(MEEP_params, PROP_params)

eps = np.abs(prop.get_epsilon())
ms = prop.meep_sim

if [False,True][1]:
    from mayavi import mlab
    dd = mlab.contour3d(eps, contours=75)
    import pdb;pdb.set_trace()

else:

    #fraction of plot lines
    ff = [0.5,0.5,0.5]

    def plot(eps, ff):
        fig, axes = plt.subplots(3,1,figsize=(11,8.5),tight_layout=True)
        #shapes
        dd = [ff[j]*eps.shape[j] for j in range(3)]
        shps = [eps.shape[j] for j in range(3)]

        #ordinate and abcissa indices
        abcs = [2,2,1]
        ords = [1,0,0]

        for j in range(3):
            axes[j].imshow(eps.take(int(dd[j]),axis=j), vmax=105, interpolation='none', \
                extent=[0,shps[abcs[j]], shps[ords[j]],0])
            axes[j].axvline(shps[abcs[j]]/2, color='k', linestyle=':')
            axes[j].axhline(shps[ords[j]]/2, color='k', linestyle=':')
            axes[j].set_xlabel(['x','y','z'][abcs[j]])
            axes[j].set_ylabel(['x','y','z'][ords[j]])
            axes[j].axvline(dd[abcs[j]], color='g')
            axes[j].axhline(dd[ords[j]], color='g')

        #Gap width
        res, wid = MEEP_params['resolution'], MEEP_params['gap_width']
        axes[2].axvline(eps.shape[1]/2-wid*res/2, color='r')
        axes[2].axvline(eps.shape[1]/2+wid*res/2, color='r')

        #Corner length
        axes[0].axvline((ms.pmlz+ms.padz)*res, color='r')
        axes[1].axvline((ms.pmlz+ms.padz)*res, color='r')

        return axes

    axes=plot(eps, ff)

    # plt.figure()
    # plt.imshow(eps[50], interpolation='none', extent=[0,eps.shape[2], eps.shape[1],0])

    import pdb;pdb.set_trace()
