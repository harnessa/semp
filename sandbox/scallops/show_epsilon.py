import numpy as np
import semp
import matplotlib.pyplot as plt;plt.ion()

MEEP_params = {

    ### Lab Properties  ###
    'waves':            0.641,

    ### Mask Properties ###
    'sim_geometry':     'edge',
    'seam_dark':        0,
    'seam_lite':        0,

    'wafer_material':   'cSi',
    'skin_material':    'metal',
    'wafer_thick':      2.,
    'skin_thick':       0.25,

    'scallop_height':   0.6,
    'scallop_depth':    0.2,
    'taper_angle':      0,
    'shave_angle':      6,

    ### Numerics ###
    'resolution':       100,
    'pml_all':          1,
    'pad_all':          0,

}

#Load propagator + analyzer
prop = semp.Propagator(MEEP_params, {'verbose':False}, is_analysis=True)
alz = semp.analysis.Analyzer({}, prop=prop, build_geo=True)

#Show epsilon
axes = alz.show_epsilon(with_lines=False)

breakpoint()
