import numpy as np
import h5py
import matplotlib.pyplot as plt;plt.ion()

# session = 'gw22_lw75_st50'
session = 'gw11_lw10_st50'
pol = ['s','p'][0]
wind = 0
is_braunbek = [False, True][1]

base_dir = f'/home/aharness/Research/Optics_Modeling/Semp_Results/full_wafer'

#Load data
fld_nme = {'s':'ez', 'p':'ey'}[pol]
with h5py.File(f'{base_dir}/{session}/fields_{pol}.h5', 'r') as f:
    fld = f[f'{fld_nme}_{wind}.r'][()] + 1j*f[f'{fld_nme}_{wind}.i'][()]

with h5py.File(f'{base_dir}/{session}/vac-fields_{pol}.h5', 'r') as f:
    vac = f[f'{fld_nme}_{wind}.r'][()] + 1j*f[f'{fld_nme}_{wind}.i'][()]

with h5py.File(f'{base_dir}/{session}/meta.h5', 'r') as f:
    xx = f['xx'][()]
    yy = f['yy'][()]
    waves = f['waves'][()]
    gap_width = f['gap_width'][()]
    device_bottom = f['device_bottom'][()]
    support_bottom = f['support_bottom'][()]

#Normalize
fld /= vac

#Subtract incident for Braunbek
if is_braunbek:
    fld -= np.heaviside(gap_width/2 + yy, 0)*np.heaviside(gap_width/2 - yy, 0)

#Extract fields at bottom of wafer
dev_fld = fld[np.argmin(abs(xx - device_bottom))]
sup_fld = fld[np.argmin(abs(xx - support_bottom))]

# plt.imshow(abs(fld), vmax=1.)
plt.imshow(np.log10(abs(fld)))

plt.figure()
plt.plot(yy, abs(dev_fld), label='Device')
plt.plot(yy, abs(sup_fld), label='Support')

plt.figure()
plt.plot(yy, np.angle(dev_fld), label='Device')
plt.plot(yy, np.angle(sup_fld), label='Support')


breakpoint()
