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

#Plot image
# axes1 = alz.show_image('hx', is_bbek=is_bbek, is_phase=is_phase)
# axes1 = alz.show_image('hy', is_bbek=is_bbek, is_phase=is_phase)
# axes1 = alz.show_image('hz', is_bbek=is_bbek, is_phase=is_phase)
# axes1 = alz.show_image('ex', is_bbek=is_bbek, is_phase=is_phase)
# axes1 = alz.show_image('ey', is_bbek=is_bbek, is_phase=is_phase)
# axes1 = alz.show_image('ez', is_bbek=is_bbek, is_phase=is_phase)

dd = abs(alz.get_data('ez'))

plt.figure()
plt.imshow(dd)

plt.figure()
plt.plot(dd[599])
plt.plot(dd[600])
plt.plot(dd[601])
plt.plot(dd[602])

axes2 = alz.show_slice('ez', is_bbek=is_bbek, is_phase=is_phase)

breakpoint()
