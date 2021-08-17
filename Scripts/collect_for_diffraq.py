"""
collect_for_diffraq.py

Author: Anthony Harness
Affiliation: Princeton University
Created on: 08-11-2021
Package: SEMP

Description: Collect and store Braunbek fields to be used by DIFFRAQ
License: Refer to $pkg_home_dir/LICENSE
"""

import numpy as np
import semp
import h5py

base_dir = 'comp_M12P3'
waves = [641, 660, 699, 725][:1]

do_save = [False, True][1]
mask = 'M12P3'
save_dir = './'
save_ext = ''

#Loop through wavelengths and collect data
data = []
for wv in waves:

    params = {
        'session':      f'{base_dir}/{wv:.0f}nm'
    }

    #Load analyzer
    alz = semp.analysis.Analyzer(params)

    #Collect and store Braunbek fields (returns sfld, pfld, xx)
    data.append(alz.collect_braunbek())

#Save data
if do_save:

    if save_ext != '':
        save_ext = '_' + save_ext

    with h5py.File(f'{save_dir}/{mask}{save_ext}.h5', 'w') as f:
        f.create_dataset('waves', data=np.array(waves)*1e-9)
        for i in range(len(waves)):
            #Write out edges
            f.create_dataset(f'{waves[i]:.0f}_s', data=data[i][0])
            f.create_dataset(f'{waves[i]:.0f}_p', data=data[i][1])
            f.create_dataset(f'{waves[i]:.0f}_x', data=data[i][2]*1e-6)
