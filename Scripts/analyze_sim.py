"""
run_sim.py

Author: Anthony Harness
Affiliation: Princeton University
Created on: 08-11-2021
Package: SEMP

Description: Run SEMP simulation
License: Refer to $pkg_home_dir/LICENSE
"""

import numpy as np
import semp
import matplotlib.pyplot as plt;plt.ion()

params = {
    # 'session':        'sommer',
    'session':          'Si_2',
    'obs_distance':     0,
}

#Load analyzer
alz = semp.analysis.Analyzer(params)

is_bbek = False

#Plot image
axes1 = alz.show_image('ez', is_bbek=is_bbek)
axes2 = alz.show_slice('ez', is_bbek=is_bbek)

breakpoint()
