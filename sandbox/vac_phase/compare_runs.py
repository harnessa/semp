import numpy as np
import matplotlib.pyplot as plt;plt.ion()
import h5py

session1 = 'vac_phase/test_edge_abs'
session2 = 'vac_phase/test_edge_pml'

pol = ['s', 'p'][1]

data_dir = '/home/aharness/Research/Optics_Modeling/Semp_Results'

comp = {'s':'ez', 'p':'hz'}[pol]

with h5py.File(f'{data_dir}/{session1}/fields_{pol}.h5', 'r') as f:
    fld1 = f[f'{comp}_0.r'][1:-1][:,1:-1] + 1j*f[f'{comp}_0.i'][1:-1][:,1:-1]

with h5py.File(f'{data_dir}/{session2}/fields_{pol}.h5', 'r') as f:
    fld2 = f[f'{comp}_0.r'][1:-1][:,1:-1] + 1j*f[f'{comp}_0.i'][1:-1][:,1:-1]

with h5py.File(f'{data_dir}/{session1}/coords_waves.h5', 'r') as f:
    xx1 = f['xx'][1:-1]
    yy1 = f['yy'][1:-1]

with h5py.File(f'{data_dir}/{session2}/coords_waves.h5', 'r') as f:
    xx2 = f['xx'][1:-1]
    yy2 = f['yy'][1:-1]

ifig, iaxes = plt.subplots(2, figsize=(8,11))
iaxes[0].imshow(abs(fld1))
iaxes[1].imshow(abs(fld2))

fig, axes = plt.subplots(2, figsize=(8,11))
axes[0].plot(xx1, abs(fld1[:,fld1.shape[1]//2]))
axes[0].plot(xx2, abs(fld2[:,fld2.shape[1]//2]), '--')

axes[1].plot(xx1, np.angle(fld1[:,fld1.shape[1]//2]))
axes[1].plot(xx2, np.angle(fld2[:,fld2.shape[1]//2]), '--')

breakpoint()
