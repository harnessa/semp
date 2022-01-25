import numpy as np
import matplotlib.pyplot as plt;plt.ion()
import h5py
import scipy.interpolate as interp
from scipy.ndimage import map_coordinates, rotate
import cv2

fnames = ['s_n0', 's_p10']

with h5py.File('s_n0.h5', 'r') as f:
    fld0 = f['field'][()]
    x0 = f['x'][()]
    y0 = f['y'][()]
    ang0 = f['rot_angle'][()]

with h5py.File('s_p10.h5', 'r') as f:
    fld1 = f['field'][()]
    ang1 = f['rot_angle'][()]

#Calculate rotation point
dpad, dpml, gap, resolution = 2, 3, 4, 10
ixc = 0.4*(dpad+dpml)*resolution
iyc = fld0.shape[0]/2
xc = 0
yc = -0.5*(dpml*2 + dpad*2 + gap) + 0.4*(dpad+dpml)

#Rotate coords
x1 = (x0 - xc)[:,None]*np.cos(ang1) - (y0 - yc)*np.sin(ang1) + xc
y1 = (x0 - xc)[:,None]*np.sin(ang1) + (y0 - yc)*np.cos(ang1) + yc

#TODO: use complex
fld0 = abs(fld0)
fld1 = abs(fld1)

#Calcluate fld1 by rotating fld0
rot_mat = cv2.getRotationMatrix2D((ixc, iyc), np.degrees(ang1), 1)
fld1 = cv2.warpAffine(fld0, rot_mat, fld0.shape[1::-1], flags=cv2.INTER_LINEAR)

plt.figure()
plt.imshow(abs(fld0))

plt.figure()
plt.imshow(abs(fld1))


thick = 1

############################################

#Select xr slice
xr = thick/2
# yr = y0
yr = np.linspace(-5,5,100)

#Rotate y into meep frame
# xm =  xr*np.cos(ang1) + yr*np.sin(ang1)
# ym = -xr*np.sin(ang1) + yr*np.cos(ang1)
xm = (xr - xc)*np.cos(ang1) - (yr - yc)*np.sin(ang1) + xc
ym = (xr - xc)*np.sin(ang1) + (yr - yc)*np.cos(ang1) + yc


#Interpolate
i1 = interp.RegularGridInterpolator((x0, y0), fld1, bounds_error=False, fill_value=0)

#Get image
# img1 = i1((x1, y1))

rot_mat = cv2.getRotationMatrix2D((ixc, iyc), -np.degrees(ang1), 1)
img1 = cv2.warpAffine(fld1, rot_mat, fld1.shape[1::-1], flags=cv2.INTER_LINEAR)


# ix1 = np.unravel_index(np.argmin(abs(x1 - xr)), x1.shape)[0]
ix1 = np.argmin(abs(x0-xr))

#Get slice
nf1 = i1((xm, ym))
# nf1 = i1((x1[ix1], y1[ix1]))

ixm = ym*resolution + fld0.shape[1]/2
iym = xm*resolution + fld0.shape[0]/2
plt.plot(ixm, iym, 'r')

nf0 = fld0[np.argmin(abs(x0-xr))]

plt.figure()
plt.plot(y0, abs(nf0))
plt.plot(ym, abs(nf1), '--')
plt.plot(y1[ix1], abs(img1[ix1]), '-.')
# plt.plot(y1[np.argmin(abs(x0-xr))], abs(img1[np.argmin(abs(x0-xr))]), '-.')


plt.figure()
plt.imshow(abs(img1))
plt.axhline(ix1, linestyle=':', color='r')


plt.figure()
plt.imshow(abs(img1 - fld0))

breakpoint()
############################################


    ############################################

xc0, yc0 = np.indices(fld0.shape)
xc1 =  xc0*np.cos(ang1) + yc0*np.sin(ang1)
yc1 = -xc0*np.sin(ang1) + yc0*np.cos(ang1)

nf2 = map_coordinates(fld1, (xc1, yc1), order=5)

plt.figure()
plt.imshow(abs(nf2))

    ############################################


breakpoint()
