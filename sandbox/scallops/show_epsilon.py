import numpy as np
import semp
import matplotlib.pyplot as plt;plt.ion()

MEEP_params = {

    ### Lab Properties  ###
    'wave':             0.641,

    ### Mask Properties ###
    'sim_geometry':     'edge',
    'seam_dark':        2,
    'seam_lite':        2,

    'wafer_material':   'cSi',
    'skin_material':    'metal',
    'wafer_thick':      2.,
    'skin_thick':       0.25,

    'scallop_height':   0.6,
    'scallop_depth':    0.1,
    'scallop_start':    0.,
    'taper_angle':      5,

    ### Numerics ###
    'resolution':       40,
    'pml_all':          1,
    'pad_all':          1,

}

#Load propagator + analyzer
prop = semp.Propagator(MEEP_params, {'verbose':False}, is_analysis=True)
alz = semp.analysis.Analyzer({}, prop=prop, build_geo=True)

#Show epsilon
axes = alz.show_epsilon(with_lines=False)

breakpoint()
