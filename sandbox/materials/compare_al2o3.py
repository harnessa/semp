import numpy as np
import matplotlib.pyplot as plt;plt.ion()
import h5py

ext1 = 'al'
ext2 = 'al2o3'

with h5py.File(f'./saves/{ext1}.h5', 'r') as f:
    yy1 = f['yy'][()]
    xind1 = f['xind'][()]
    ez1 = f['ez'][xind1]
    hz1 = f['hz'][xind1]

with h5py.File(f'./saves/{ext2}.h5', 'r') as f:
    yy2 = f['yy'][()]
    xind2 = f['xind'][()]
    ez2 = f['ez'][xind2]
    hz2 = f['hz'][xind2]

#Plot
fig, axes = plt.subplots(2, 2, figsize=(9,9))
axes[0,0].plot(yy1, abs(ez1), '-',  label='Al')
axes[0,0].plot(yy2, abs(ez2), '--', label='Al203')

axes[1,0].plot(yy1, abs(hz1), '-')
axes[1,0].plot(yy2, abs(hz2), '--')


axes[0,1].plot(yy1, np.angle(ez1), '-')
axes[0,1].plot(yy2, np.angle(ez2), '--')

axes[1,1].plot(yy1, np.angle(hz1), '-')
axes[1,1].plot(yy2, np.angle(hz2), '--')

axes[0,0].legend()

breakpoint()
