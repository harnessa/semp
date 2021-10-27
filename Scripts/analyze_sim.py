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
    'session':          'final_model/m12px_gaps_w223',
    'obs_distance':     0.,
}

#Load analyzer
alz = semp.analysis.Analyzer(params)

is_bbek = True
is_phase = False
wave = 0.641

#Plot image
img, axes1 = alz.show_image('hz', wave=wave, is_bbek=is_bbek, is_phase=is_phase)#, vmax=[1.25,None][int(is_phase)])
# img, axes1 = alz.show_image('ez', wave=wave, is_bbek=is_bbek, is_phase=True)#, vmax=1.25)
slc, axes2 = alz.show_slice('hz', wave=wave, is_bbek=is_bbek, is_phase=is_phase)

# axes1.set_xlim([-4*wave,4*wave])
# axes1.set_ylim([ 4*wave,-4*wave])

breakpoint()
