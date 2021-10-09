import numpy as np
import h5py
import matplotlib.pyplot as plt;plt.ion()


with h5py.File('./normal_ez.h5', 'r') as f:
    nez = f['ez'][()]

with h5py.File('./dft_ez.h5', 'r') as f:
    dez = f['ez'][()]

plt.imshow(abs(abs(nez)/abs(nez).max() - abs(dez)/abs(dez).max()))
print(abs(abs(nez)/abs(nez).max() - abs(dez)/abs(dez).max()).max())

breakpoint()
