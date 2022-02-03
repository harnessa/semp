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
    'session':          'skin_epsi_1',
    # 'session':          'test_all2',
    # 'base_dir':         f'{semp.tmp_dir}/tests',
    # 'session':          'thick_screen_all',
    'obs_distance':     0.,
}

#Load analyzer
alz = semp.analysis.Analyzer(params)

is_bbek = True
is_phase = False
wave = 0.641

plt.ioff()
#Plot image
# img, axes1 = alz.show_image('ez', wave=wave, is_bbek=is_bbek, is_phase=is_phase)#, vmax=[1.25,None][int(is_phase)])
# img, axes1 = alz.show_image('ez', wave=wave, is_bbek=is_bbek, is_phase=True)#, vmax=1.25)
slc1, axes2 = alz.show_slice('ez', wave=wave, is_bbek=is_bbek, is_phase=is_phase)
slc3, axes2 = alz.show_slice('ez', wave=wave, is_bbek=is_bbek, is_phase=not is_phase)

params['session'] = 'skin_epsi_2'
alz2 = semp.analysis.Analyzer(params)

slc2, axes3 = alz2.show_slice('ez', wave=wave, is_bbek=is_bbek, is_phase=is_phase)
slc4, axes3 = alz2.show_slice('ez', wave=wave, is_bbek=is_bbek, is_phase=not is_phase)

plt.ion()

plt.figure()
plt.plot(slc1)
plt.plot(slc2, '--')

plt.figure()
plt.plot(slc3)
plt.plot(slc4, '--')

breakpoint()

# axes1.set_xlim([-4*wave,4*wave])
# axes1.set_ylim([ 4*wave,-4*wave])



# yy1 = alz.yy.copy()
#
# params['session'] = 'final_model/22gap'
# alz = semp.analysis.Analyzer(params)
# slc2, axes3 = alz.show_slice('ez', wave=wave, is_bbek=is_bbek, is_phase=is_phase)
# yy2 = alz.yy.copy()
#
#
# plt.figure()
# plt.plot(yy1, slc**2)
# plt.plot(yy2, slc2**2)
#
# breakpoint()
