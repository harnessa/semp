import numpy as np
import h5py
import matplotlib.pyplot as plt;plt.ion()


with h5py.File('./saves/si_broad_ez_noinjHy.h5', 'r') as f:
    waves = f['waves'][()]
    nez = f['ez'][()]
    nhy = f['hy'][()]

with h5py.File('./saves/si_broad_ez_injHy.h5', 'r') as f:
    bez = f['ez'][()]
    bhy = f['hy'][()]

fig, axes = plt.subplots(1, 3, figsize=(11,6))
for i in range(len(waves)):
    cn = abs(nez[i])#/abs(nez[i]).max()
    cb = abs(bez[i])#/abs(bez[i]).max()
    axes[0].imshow(cn)
    axes[1].imshow(cb)
    # axes[2].imshow(np.log10(abs(cn-cb)))
    print(abs(nez[i] - bez[i]).max())
    print(abs(nhy[i] - bhy[i]).max())

    breakpoint()
breakpoint()
