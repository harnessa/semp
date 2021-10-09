import numpy as np
import h5py
import matplotlib.pyplot as plt;plt.ion()

with h5py.File('./edge_ez.h5', 'r') as f:
    fez = f['ez'][()]


with h5py.File('./edge_time_ez.h5', 'r') as f:
    tez = f['ez'][()]


plt.imshow(abs(fez))

plt.figure()
plt.imshow(abs(abs(fez)/abs(fez).max() - abs(tez)/abs(tez).max()), vmax=0.01)

print(abs(abs(fez)/abs(fez).max() - abs(tez)/abs(tez).max())[100:].max())

breakpoint()
