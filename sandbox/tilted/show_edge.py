import numpy as np
import matplotlib.pyplot as plt;plt.ion()
import h5py

fnames = ['s_n0', 's_n10', 's_p10'][:]

data, xxs, yys, angs = [], [], [], []
for fn in fnames:
    with h5py.File(fn + '.h5', 'r') as f:
        data.append(f['field'][()])
        x = f['x'][()]
        y = f['y'][()]
        ang = f['rot_angle'][()]

    #Rotate coords
    nx =  x[:,None]*np.cos(ang) + y*np.sin(ang)
    ny = -x[:,None]*np.sin(ang) + y*np.cos(ang)

    xxs.append(nx)
    yys.append(ny)
    angs.append(ang)

for i in range(len(fnames)):
    plt.figure()
    plt.imshow(abs(data[i]))

breakpoint()
# #Interpolate rotation
# from scipy.interpolate import RectBivariateSpline
# x2, y2 = np.arange(len(x)), np.arange(len(y))
# newd = []
# for i in range(len(data)):
#     interpolator = RectBivariateSpline(x,y, abs(data[i]))
#
#     breakpoint()

# breakpoint()

#Derotate
from scipy.ndimage import rotate
newd2, y2 = [], []
for i in range(len(data)):
    newd2.append(rotate(data[i], np.degrees(-angs[i]), reshape=False, order=5))
    y2.append(rotate(yys[i], np.degrees(-angs[i]), reshape=False, order=5))

from skimage.transform import AffineTransform, warp
newd3, y3 = [], []
for i in range(len(data)):
    tran = AffineTransform(rotation=-angs[i])
    fr = warp(data[i].real, tran)
    fi = warp(data[i].imag, tran)
    newd3.append(fr + 1j*fi)
    y3.append(warp(yys[i], tran))

for i in range(len(newd3)):
    plt.figure()
    plt.imshow(abs(newd3[i]))

thick = 1
xind = np.where(xxs[0][:,0] > thick/2)[0][0]
plt.figure()
plt.plot(y3[0][xind], np.angle(data[0][xind]), label=f'{np.degrees(angs[0]):.0f}')
plt.plot(y3[1][xind], np.angle(newd3[1][xind]), '--', label=f'{np.degrees(angs[1]):.0f}')
plt.plot(y3[2][xind], np.angle(newd3[2][xind]), '-.', label=f'{np.degrees(angs[2]):.0f}')
plt.legend()

# plt.figure()
# plt.plot(abs(data[0][xind]))
# plt.plot(abs(newd2[1][xind]), '--')
# plt.plot(abs(newd3[1][xind]), '-.')
#
# plt.figure()
# plt.plot(abs(data[0][xind]))
# plt.plot(abs(newd2[2][xind]), '--')
# plt.plot(abs(newd3[2][xind]), '-.')

breakpoint()
