import numpy as np
import h5py
import cv2
import matplotlib.pyplot as plt;plt.ion()

def get_slice(ang, pol, wind, base_dir):

    sgn = ["n", "p"][int((np.sign(ang)+1)/2)]
    load_dir = f'{base_dir}/{sgn}{abs(ang)*100:.0f}'

    #Load meta data
    with h5py.File(f'{load_dir}/meta.h5', 'r') as f:
        xx = f['xx'][()]        #Source coordinates (light source)
        yy = f['yy'][()]
        waves = f['waves'][()]
        rot_angle = f['rot_angle'][()]
        wafer_center = f['wafer_center'][()]
        # wafer_thick = f['wafer_thick'][()]
        # resolution = f['resolution'][()]
        # edge_y = f['edge_y'][()]

    #FIXME:
    wafer_thick = 2.25
    resolution = 30
    edge_y = -12.3

    #Load data
    comp = {'s':'ez', 'p':'ey'}[pol]
    with h5py.File(f'{load_dir}/fields_{pol}.h5', 'r') as f:
        waf_fld = f[f'{comp}_{wind}.r'][()] + 1j*f[f'{comp}_{wind}.i'][()]

    with h5py.File(f'{load_dir}/vac-fields_{pol}.h5', 'r') as f:
        vac_fld = f[f'{comp}_{wind}.r'][()] + 1j*f[f'{comp}_{wind}.i'][()]

    #Normalize
    waf_fld /= vac_fld

    #Turn into Braunbek field
    waf_fld -= np.heaviside(yy - edge_y, 0)

    #Get rotation center in pixel coordinates
    rot_cen = (wafer_center * resolution + np.array(waf_fld.shape)/2)[::-1]

    #Derotate image
    rot_mat = cv2.getRotationMatrix2D(tuple(rot_cen), -np.degrees(rot_angle), 1)
    do_rot = lambda fld:  \
        cv2.warpAffine(fld, rot_mat, waf_fld.shape[1::-1], flags=cv2.INTER_LINEAR)

    newr = do_rot(waf_fld.real)
    newi = do_rot(waf_fld.imag)
    new_fld = newr + 1j*newi

    #Get slice at bottom of wafer
    xind = np.argmin(abs(xx - wafer_thick/2))
    slc = new_fld[xind]
    sy = yy - edge_y


    # plt.figure()
    # plt.imshow(abs(waf_fld))
    #
    # plt.figure()
    # plt.imshow(abs(new_fld))
    #
    # fig, axes = plt.subplots(2, sharex=True, figsize=(8,11))
    # axes[0].plot(sy, abs(slc))
    # axes[1].plot(sy, np.angle(slc))
    #
    # breakpoint()

    return sy, slc

############################################
####	Get Slices ####
############################################

#Angles to run
ang_max = 5
nangs = 11
session = 'tilted_runs'

base_dir = f'/home/aharness/Research/Optics_Modeling/Semp_Results/{session}'
waves = [641, 660, 699, 725]

#Angles
nangs += (nangs + 1) % 2
angs = np.linspace(-ang_max, ang_max, nangs)

#Get slices
sdata, pdata, ydata = [], [], []
for wind in range(len(waves)):
    stmp, ptmp, ytmp = [], [], []
    for ang in angs:

        #Get data
        sy, sfld = get_slice(ang, 's', wind, base_dir)
        py, pfld = get_slice(ang, 'p', wind, base_dir)

        #Append
        stmp.append(sfld)
        ptmp.append(pfld)
        ytmp.append(sy)

    #Append
    sdata.append(np.array(stmp))
    pdata.append(np.array(ptmp))
    ydata.append(np.array(ytmp))

#Save data
if True:
    with h5py.File(f'./Results/results__{session}.h5', 'w') as f:
        f.create_dataset('waves', data=np.array(waves)*1e-9)
        f.create_dataset('angles', data=np.array(angs))
        for i in range(len(waves)):
            #Write out edges
            f.create_dataset(f'{waves[i]:.0f}_x', data=ydata[i]*1e-6)
            f.create_dataset(f'{waves[i]:.0f}_s', data=sdata[i])
            f.create_dataset(f'{waves[i]:.0f}_p', data=pdata[i])
