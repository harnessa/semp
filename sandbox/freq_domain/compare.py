import numpy as np
import h5py
import matplotlib.pyplot as plt;plt.ion()


with h5py.File('./saves/si_narrow_ez.h5', 'r') as f:
    waves = f['waves'][()]
    nez = f['ez'][()]

with h5py.File('./saves/si_solver.h5', 'r') as f:
    bez = f['ez'][()]

fig, axes = plt.subplots(1, 3, figsize=(11,6))
for i in range(len(waves)):
    cn = abs(nez[i])/abs(nez[i]).max()
    cb = abs(bez[i])/abs(bez[i]).max()
    axes[0].imshow(cn)
    axes[1].imshow(cb)
    axes[2].imshow(abs(cn-cb))
    print(abs(cn - cb).max())

    breakpoint()
breakpoint()
