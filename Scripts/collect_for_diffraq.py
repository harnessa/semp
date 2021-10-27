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
import matplotlib.pyplot as plt;plt.ion()

session = 'final_model/m12px'
waves = [641, 660, 699, 725]

do_save = [False, True][1]
mask = 'M12PX'
save_ext = 'final_test'
save_dir = '/home/aharness/repos/diffraq/External_Data/Vector_Edges'

#Figures
fig, axes = plt.subplots(2, sharex=True, sharey=True, figsize=(6,9))


#Load analyzer
params = {
    'session':      session,
}
alz = semp.analysis.Analyzer(params)

#Loop through wavelengths and collect data
data = []
for wv in waves:

    #Collect and store Braunbek fields (returns sfld, pfld, xx)
    data.append(alz.collect_braunbek(wave=wv*1e-3))

    #Plot
    axes[0].plot(data[-1][2], np.abs(data[-1][0]), label=f'{wv}nm')
    axes[1].plot(data[-1][2], np.abs(data[-1][1]), label=f'{wv}nm')

axes[0].legend()

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

else:
    breakpoint()
