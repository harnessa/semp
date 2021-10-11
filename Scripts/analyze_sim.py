"""
run_sim.py

Author: Anthony Harness
Affiliation: Princeton University
Created on: 08-11-2021
Package: SEMP

Description: Analyze SEMP simulation
License: Refer to $pkg_home_dir/LICENSE
"""

import numpy as np
import semp
import matplotlib.pyplot as plt;plt.ion()

params = {
    'session':          'M12P6_h6_d1_a1/641nm',
    'obs_distance':     0,
}

#Load analyzer
alz = semp.analysis.Analyzer(params)

is_bbek = False
is_phase = False
wave = 0.641

#Plot image
axes1 = alz.show_image('ez', wave=wave, is_bbek=is_bbek, is_phase=is_phase)
axes2 = alz.show_slice('ez', wave=wave, is_bbek=is_bbek, is_phase=is_phase)

breakpoint()
