import numpy as np
import h5py
import matplotlib.pyplot as plt;plt.ion()


with h5py.File('./saves/si_broad_ez_pml8.h5', 'r') as f:
    waves = f['waves'][()]
    nez = f['ez'][()]

with h5py.File('./saves/si_broad_ez_pml12.h5', 'r') as f:
    bez = f['ez'][()]

fig, axes = plt.subplots(1, 3, figsize=(11,6))
for i in range(len(waves)):
    cn = abs(nez[i])/abs(nez[i]).max()
    cb = abs(bez[i])/abs(bez[i]).max()
    dd = cb.shape[0]-cn.shape[0]
    cb = cb[dd//2:-dd//2][:,dd//2:-dd//2]
    axes[0].imshow(cn)
    axes[1].imshow(cb)
    axes[2].imshow(np.log10(abs(cn-cb)))
    print(abs(cn - cb).max())

    breakpoint()
breakpoint()
