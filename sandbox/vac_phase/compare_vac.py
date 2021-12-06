import numpy as np
import matplotlib.pyplot as plt;plt.ion()
import h5py

session = 'vac_phase/test_pml'

pol = ['s', 'p'][0]

data_dir = '/home/aharness/Research/Optics_Modeling/Semp_Results'

comp = {'s':'ez', 'p':'hz'}[pol]

with h5py.File(f'{data_dir}/{session}/fields_{pol}.h5', 'r') as f:
    fld = f[f'{comp}_0.r'][1:-1][:,1:-1] + 1j*f[f'{comp}_0.i'][1:-1][:,1:-1]

with h5py.File(f'{data_dir}/{session}/vac-fields_{pol}.h5', 'r') as f:
    vac = f[f'{comp}_0.r'][1:-1] + 1j*f[f'{comp}_0.i'][1:-1]

with h5py.File(f'{data_dir}/{session}/coords_waves.h5', 'r') as f:
    xx = f['xx'][1:-1]
    yy = f['yy'][1:-1]

plt.figure()
plt.imshow(abs(fld))

fig, axes = plt.subplots(2, figsize=(8,11), sharex=True)
axes[0].plot(xx, abs(fld[:,fld.shape[1]//2]))
axes[0].plot(xx, abs(vac), '--')

axes[1].plot(xx, np.angle(fld[:,fld.shape[1]//2]))
axes[1].plot(xx, np.angle(vac), '--')

breakpoint()
