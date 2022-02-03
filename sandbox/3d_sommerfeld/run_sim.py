import numpy as np
import meep as mp
import h5py
from scipy.special import fresnel
import cv2
import matplotlib.pyplot as plt;plt.ion()

############################################
####	User Input ####
############################################

resolution = 30

pol = ['s','p'][0]

inn_rot_angle = np.radians(0)
out_rot_angle = np.radians(0)

dpml = 2
dpad = 4
gap = 6
thick = 0.05#min(0.01, 0.5/resolution)

wave = 0.6
df = 0.02  # turn-on bandwidth

############################################
############################################

#Derived
# fcen = 1/wave
fcen = 1/wave * np.cos(out_rot_angle)
lx_np = dpad*2 + thick
ly_np = dpad*2 + gap
cell_size = mp.Vector3(dpml*2 + lx_np, dpml*2 + ly_np)
src_comp = getattr(mp, {'s': 'Ez', 'p': 'Ey'}[pol])
kk = 2*np.pi/wave

#Kpoint
# kpoint = mp.Vector3(z=np.sin(out_rot_angle)).scale(fcen)
kpoint = mp.Vector3()

kz_2d = ['complex', 'real/imag'][0]

def pw_amp(k, x0):
    def _pw_amp(x):
        return np.exp(1j * k.dot(x + x0))
    return _pw_amp

#Source
sim_src = mp.GaussianSource(fcen, fwidth=df, is_integrated=True)
sources = [mp.Source(sim_src, component=src_comp,
            center=mp.Vector3(x=-0.5*lx_np), \
            size=mp.Vector3(y=cell_size.y), )]
            # amp_func=pw_amp(kpoint, mp.Vector3(x=-0.5*lx_np)))]

#PML
if pol == 's':
    pml_layers = [mp.PML(thickness=dpml)]
else:
    pml_layers = [mp.PML(thickness=dpml,  direction=mp.X), \
        mp.PML(thickness=dpml, direction=mp.Y, side=-1),
        mp.Absorber(thickness=dpml, direction=mp.Y, side=1)]

#Geometry
geometry = []
waf_sy = dpad+dpml
e1 = mp.Vector3(x=1).rotate(mp.Vector3(z=1), inn_rot_angle)
e2 = mp.Vector3(y=1).rotate(mp.Vector3(z=1), inn_rot_angle)
wafer = mp.Block(material=mp.metal, size=mp.Vector3(thick, waf_sy, mp.inf),
    center=mp.Vector3(y=-0.5*cell_size.y + 0.4*waf_sy), e1=e1, e2=e2)
geometry += [wafer]

edge_y = wafer.center.y + wafer.size.y/2

############################################
############################################

from mpi4py import MPI
size = MPI.COMM_WORLD.size

sgn = ["n", "p"][int((np.sign(inn_rot_angle)+1)/2)]
ext = f'{sgn}{abs(np.degrees(inn_rot_angle)):.0f}_{abs(np.degrees(out_rot_angle)):.0f}'

if size > 1:

    #Build simulation
    sim = mp.Simulation(force_complex_fields=True,
        resolution=resolution, ensure_periodicity=False,
        cell_size=cell_size, boundary_layers=[mp.PML(thickness=dpml)], sources=sources,
        geometry=[], k_point=kpoint, kz_2d=kz_2d)

    #Add DFT fields
    vol = mp.Volume(size=mp.Vector3(lx_np, ly_np))
    dft_obj = sim.add_dft_fields([src_comp], [fcen], where=vol)

    #Run vacuum
    sim.run(until_after_sources=mp.stop_when_fields_decayed(50, src_comp, mp.Vector3(-0.5*cell_size.x+dpml), 1e-6))
    vac = sim.get_dft_array(dft_obj, src_comp, 0)

    #Reset with geometry
    sim.reset_meep()
    sim = mp.Simulation(force_complex_fields=True,
        resolution=resolution, ensure_periodicity=False,
        cell_size=cell_size, boundary_layers=pml_layers, sources=sources,
        geometry=geometry, k_point=kpoint, kz_2d=kz_2d)
    dft_obj = sim.add_dft_fields([src_comp], [fcen], where=vol)

    #Run
    sim.run(until_after_sources=mp.stop_when_fields_decayed(50, src_comp, mp.Vector3(-0.5*cell_size.x+dpml), 1e-6))

    #Get fields
    fld = sim.get_dft_array(dft_obj, src_comp, 0)
    x,y,z,w = sim.get_array_metadata(dft_cell=dft_obj)

    if mp.am_master():
        with h5py.File(f'./saved_{ext}.h5', 'w') as f:
            f.create_dataset('fld', data=fld)
            f.create_dataset('vac', data=vac)
            f.create_dataset('x', data=x)
            f.create_dataset('y', data=y)
    import sys;sys.exit()

else:

    with h5py.File(f'./saved_{ext}.h5', 'r') as f:
        fld = f['fld'][()]
        vac = f['vac'][()]
        x = f['x'][()]
        y = f['y'][()]

#Normalize
fld /= vac


# ext1 = 'test_k'
# ext2 = 'test_kb'
#
# with h5py.File(f'./saved_{ext1}.h5', 'r') as f:
#     fld1 = f['fld'][()]
#     vac1 = f['vac'][()]
#     x = f['x'][()]
#     y = f['y'][()]
#
# #Normalize
# fld1 /= vac1
#
# with h5py.File(f'./saved_{ext2}.h5', 'r') as f:
#     fld2 = f['fld'][()]
#     vac2 = f['vac'][()]
#
# #Normalize
# fld2 /= vac2
# plt.imshow(abs(fld1-fld2))
# breakpoint()

############################################
####	Sommerfeld ####
############################################

def G_func(s):
    S,C = fresnel(np.sqrt(2./np.pi)*s)
    ans = (1. + 1j)/2 - (C + 1j*S)
    ans *= np.exp(-1j*s**2.)
    del S, C
    return ans

def get_sommerfeld_solution(xx, yy, phi0, beta, kk, is_bbek=False):

    #Get coordinates (x=y, y=z)
    rho = np.sqrt(yy**2 + xx**2.)

    phi = np.arctan2(xx,yy) + np.pi
    zind = np.isclose(xx,0)
    phi[zind] = (2 - np.heaviside(yy[zind], 1))*np.pi

    zz = 27e6*np.sin(beta)

    #Incident field
    I0 = np.exp(1j*kk*(yy*np.cos(phi0)*np.cos(beta) + xx*np.sin(phi0)*np.cos(beta) \
        - zz*np.sin(beta)))

    #Build arguments of calculation
    uu = -np.sqrt(2.*kk*rho*np.cos(beta))*np.cos(0.5*(phi - phi0))
    vv = -np.sqrt(2.*kk*rho*np.cos(beta))*np.cos(0.5*(phi + phi0))
    pre = (np.exp(-1j*np.pi/4.) / np.sqrt(np.pi) * np.sqrt(np.pi/2)) * \
        np.exp(1j*kk*rho*np.cos(beta)) * np.cos(beta) * np.exp(-1j*kk*zz*np.sin(beta))

    #Normalize prefactor by plane wave
    pre /= I0

    #Intermediate calculations
    Umid = pre*G_func(uu)
    Gv = pre*G_func(vv)
    Dmid = 2j*pre / np.sqrt(np.pi*kk*rho*np.cos(beta))

    #Subtract incident field if Braunbek
    if is_bbek:
        Umid -= np.heaviside(yy, 1)

    #Cleanup
    del uu, vv, pre

    #Get field solution for s,p polarization
    Ez = Umid - Gv                  #s-pol
    Hz = Umid + Gv                  #p-pol

    #Cleanup
    del Umid, Gv, Dmid, rho, phi

    return Ez, Hz, I0

phi0 = np.pi/2 + inn_rot_angle
xx = np.tile(x - thick/2, (len(y),1)).T
yy = np.tile(y - edge_y, (len(x),1))

func, vmax = [[np.abs, 0.1], [np.angle,0.5]][1]

ez, hz, i0 = get_sommerfeld_solution(xx, yy, phi0, out_rot_angle, kk)

ez0, hz, i0b = get_sommerfeld_solution(xx, yy, phi0, 0, kk*np.cos(out_rot_angle))
# diff = ez/ez0
# # fld *= diff
# diff = ez - ez0
# fld += diff

# plt.figure()
# plt.imshow(func(ez))
# plt.figure()
# plt.imshow(func(ez0))
# plt.figure()
# plt.imshow(abs(func(ez)-func(ez0)))
# # plt.imshow((func(ez)-func(ez0)))

# plt.figure()
# plt.plot(abs(ez[81]))
# plt.plot(abs(ez0[81]),'--')
#
# plt.figure()
# plt.plot(np.angle(ez[81]))
# plt.plot(np.angle(ez0[81]),'--')
# breakpoint()

############################################
############################################

#Get rotation center in pixel coordinates
if not np.isclose(inn_rot_angle, 0):
    rot_cen = (np.array([wafer.center.x, wafer.center.y]) * resolution + np.array(fld.shape)/2)[::-1]

    #Derotate image
    rot_mat = cv2.getRotationMatrix2D(tuple(rot_cen), -np.degrees(inn_rot_angle), 1)
    do_rot = lambda fld:  \
        cv2.warpAffine(fld, rot_mat, fld.shape[1::-1], flags=cv2.INTER_LINEAR)

    newr = do_rot(fld.real)
    newi = do_rot(fld.imag)
    fld = newr + 1j*newi


# plt.figure()
# plt.imshow(func(fld))
# plt.figure()
# plt.imshow(func(ez))
# plt.figure()
# plt.imshow(abs(func(fld)-func(ez))/func(ez)*100, vmax=vmax*100)
# plt.imshow(abs(func(fld)-func(ez)), vmax=vmax)
plt.imshow(abs(fld-ez), vmax=0.1)

# plt.figure()
# plt.plot(np.angle(i0)[:,200])
# plt.plot(np.angle(vac)[:,200],'--')
breakpoint()
