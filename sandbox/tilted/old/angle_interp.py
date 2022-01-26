import numpy as np
import h5py
import matplotlib.pyplot as plt;plt.ion()
from scipy.interpolate import RectBivariateSpline, interp2d, RegularGridInterpolator

session = 'tilted_runs'

wave = 641

with h5py.File(f'../Results/results__{session}.h5', 'r') as f:

    angles = f['angles'][()]
    xx = f[f'{wave}_x'][0]
    sfld = f[f'{wave}_s'][()]
    pfld = f[f'{wave}_p'][()]

xx = xx[:1200]
sfld = sfld[:,:1200]

for i in range(len(angles)):
    plt.plot(xx, abs(sfld[i]), label=f'{angles[i]:.2f}')
plt.legend()
breakpoint()

i0, i1 = 2, 3
test_ang = (angles[i0] + angles[i1])/2
# test_xx = np.linspace(-5, 20, 1000) * 1e-6
test_xx = np.linspace(xx.min(),xx.max(),500)

sintr = RectBivariateSpline(angles, xx, sfld.real, kx=5, ky=5)
sinti = RectBivariateSpline(angles, xx, sfld.imag, kx=5, ky=5)

sint2 = RegularGridInterpolator((angles, xx), sfld, bounds_error=False, fill_value=0)

sint3a = RectBivariateSpline(angles, xx, abs(sfld), kx=5, ky=5)
sint3p = RectBivariateSpline(angles, xx, np.angle(sfld), kx=5, ky=5)
f3 = (sint3a(test_ang, test_xx) * np.exp(1j*sint3p(test_ang, test_xx)))[0]


f1 = (sintr(test_ang, test_xx) + 1j*sinti(test_ang, test_xx))[0]
f2 = sint2((test_ang, test_xx))

plt.figure()
plt.plot(xx, abs(sfld[i0]), 'k')
plt.plot(xx, abs(sfld[i1]), 'r')
plt.plot(test_xx, abs(f1), '--')
plt.plot(test_xx, abs(f2), '-.')
plt.plot(test_xx, abs(f3), ':')

plt.figure()
plt.plot(xx, np.angle(sfld[i0]), 'k')
plt.plot(xx, np.angle(sfld[i1]), 'r')
plt.plot(test_xx, np.angle(f1), '--')
plt.plot(test_xx, np.angle(f2), '-.')
plt.plot(test_xx, np.angle(f3), ':')

# plt.figure()
# plt.plot(xx, np.unwrap(np.angle(sfld[i0])), 'k')
# plt.plot(xx, np.unwrap(np.angle(sfld[i1])), 'r')
# plt.plot(test_xx, np.unwrap(np.angle(f1)), '--')
# plt.plot(test_xx, np.unwrap(np.angle(f2)), '-.')
# plt.plot(test_xx, np.unwrap(np.angle(f3)), ':')

breakpoint()
