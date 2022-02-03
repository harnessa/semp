import numpy as np
import h5py
import matplotlib.pyplot as plt;plt.ion()

# session = 'gw22_lw75_st50'
session = 'gw11_lw10_st50'

is_dev = [False, True][0]


base_dir = f'/home/aharness/Research/Optics_Modeling/Semp_Results/full_wafer'

def get_field(pol, wind, is_dev=False):

    #Load data
    fld_nme = {'s':'ez', 'p':'ey'}[pol]
    with h5py.File(f'{base_dir}/{session}/fields_{pol}.h5', 'r') as f:
        fld = f[f'{fld_nme}_{wind}.r'][()] + 1j*f[f'{fld_nme}_{wind}.i'][()]

    with h5py.File(f'{base_dir}/{session}/vac-fields_{pol}.h5', 'r') as f:
        vac = f[f'{fld_nme}_{wind}.r'][()] + 1j*f[f'{fld_nme}_{wind}.i'][()]

    with h5py.File(f'{base_dir}/{session}/meta.h5', 'r') as f:
        xx = f['xx'][()]
        yy = f['yy'][()]
        gap_width = f['gap_width'][()]
        device_bottom = f['device_bottom'][()]
        support_bottom = f['support_bottom'][()]

    #Normalize
    fld /= vac

    #Extract fields at bottom of wafer
    if is_dev:
        fld = fld[np.argmin(abs(xx - device_bottom))]
    else:
        fld = fld[np.argmin(abs(xx - support_bottom))]

    #Subtract incident for Braunbek
    fld -= np.heaviside(gap_width/2 + yy, 0)*np.heaviside(gap_width/2 - yy, 0)

    #Shift y to edge
    yy += gap_width/2

    #Cut in half
    fld = fld[:len(fld)//2]
    yy = yy[:len(yy)//2]

    return fld, yy

waves = [641, 660, 699, 725]

fig, axes = plt.subplots(1, 2)

sdata, pdata, ydata = [], [], []
for i in range(len(waves)):

    sfld, yy = get_field('s', i, is_dev=is_dev)
    pfld, yy = get_field('p', i, is_dev=is_dev)

    axes[0].plot(yy, abs(sfld), label=waves[i])
    axes[1].plot(yy, abs(pfld), label=waves[i])

    sdata.append(sfld)
    pdata.append(pfld)
    ydata.append(yy)

axes[0].legend()

if [False, True][1]:

    save_dir = '/home/aharness/repos/Milestone_2/diffraq_analysis/modeling/Vector_Edges/full_wafer/saves'
    ext = ['support', 'device'][int(is_dev)]


    with h5py.File(f'{save_dir}/{session}__{ext}.h5', 'w') as f:
        for i in range(len(waves)):
            f.create_dataset(f'{waves[i]}_s', data=sdata[i])
            f.create_dataset(f'{waves[i]}_p', data=pdata[i])
            f.create_dataset(f'{waves[i]}_x', data=ydata[i]*1e-6)


breakpoint()
