import numpy as np
import h5py
# import cv2
import matplotlib.pyplot as plt;plt.ion()

def get_data(ext, pol, wind):

    base_dir = '/home/aharness/Research/Optics_Modeling/Semp_Results/tilted_runs'
    load_dir = f'{base_dir}/{ext}'

    #Load meta data
    with h5py.File(f'{load_dir}/meta.h5', 'r') as f:
        xx = f['xx'][()]        #Source coordinates (light source)
        yy = f['yy'][()]
        waves = f['waves'][()]
        inn_rot_angle = f['inn_rot_angle'][()]
        wafer_center = f['wafer_center'][()]
        wafer_thick = f['wafer_thick'][()]
        resolution = f['resolution'][()]
        edge_y = f['edge_y'][()]

    #Load data
    comp = {'s':'ez', 'p':'ey'}[pol]
    with h5py.File(f'{load_dir}/fields_{pol}.h5', 'r') as f:
        waf_fld = f[f'{comp}_{wind}.r'][()] + 1j*f[f'{comp}_{wind}.i'][()]

    with h5py.File(f'{load_dir}/vac-fields_{pol}.h5', 'r') as f:
        vac_fld = f[f'{comp}_{wind}.r'][()] + 1j*f[f'{comp}_{wind}.i'][()]

    #Turn into Braunbek field
    # waf_fld -= np.heaviside(yy - edge_y, 0)

    # #Get rotation center in pixel coordinates
    # rot_cen = (wafer_center * resolution + np.array(waf_fld.shape)/2)[::-1]
    #
    # #Derotate image
    # rot_mat = cv2.getRotationMatrix2D(tuple(rot_cen), -np.degrees(inn_rot_angle), 1)
    # do_rot = lambda fld:  \
    #     cv2.warpAffine(fld, rot_mat, waf_fld.shape[1::-1], flags=cv2.INTER_LINEAR)
    #
    # newr = do_rot(waf_fld.real)
    # newi = do_rot(waf_fld.imag)
    # new_fld = newr + 1j*newi
    #
    # #Get slice at bottom of wafer
    # xind = np.argmin(abs(xx - wafer_thick/2))
    # slc = new_fld[xind]
    # sy = yy - edge_y
    #
    # return sy, slc

    return waf_fld, vac_fld

pol = ['s','p'][1]
wind = 0

ext = '3D_test/n0_1000'

waf_fld, vac_fld = get_data(ext, pol, wind)

#Normalize
waf_fld /= vac_fld

plt.figure()
plt.imshow(abs(waf_fld))
plt.figure()
plt.imshow(abs(vac_fld))

breakpoint()
