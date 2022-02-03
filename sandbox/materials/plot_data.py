import numpy as np
import matplotlib.pyplot as plt;plt.ion()
import meep as mp
from meep import materials as mat_lib

w0, w1 = 0.6, 0.8

def get_data(ext):
    data = np.genfromtxt(ext + '.csv', delimiter=',')

    #split n, k
    i1 = np.where(np.isnan(data[:,0]))[0][-1]
    wave = data[1:i1][:,0]
    nn = data[1:i1][:,1]
    kk = data[i1+1:][:,1]

    #permittivity
    perm = (nn**2 - kk**2) + 1j*(2*nn*kk)

    #trim
    perm = perm[(wave >= w0) & (wave <= w1)]
    wave = wave[(wave >= w0) & (wave <= w1)]

    return wave, perm

ext1 = 'Si_Aspnes'
ext2 = 'Si_Schinke'

#Get data
wave1, pexp1 = get_data(ext1)
wave2, pexp2 = get_data(ext2)

#meep data
wavem = np.linspace(w0, w1, 500)
med1 = mat_lib.cSi
perm1 = np.array([med1.epsilon(1/w)[0,0] for w in wavem])

#defined material
ceps = 14.979309348016363+0.1177371364850722j
Dcon = 2*np.pi/0.641 * ceps.imag / ceps.real

med2 = mp.Medium(epsilon=ceps.real, D_conductivity=Dcon)
perm2 = np.array([med2.epsilon(1/w)[0,0] for w in wavem])

um_scale = 1.0
# conversion factor for eV to 1/um [=1/hc]
eV_um_scale = um_scale/1.23984193
cSi_range = mp.FreqRange(min=um_scale, max=um_scale/0.4)

cSi_frq1 = 3.64/um_scale
cSi_gam1 = 0
cSi_sig1 = 8
cSi_frq2 = 2.76/um_scale
cSi_gam2 = 2*0.063/um_scale
cSi_sig2 = 2.85
cSi_frq3 = 1.73/um_scale
cSi_gam3 = 2*2.5/um_scale
cSi_sig3 = -0.107

cSi_susc = []
cSi_susc += [mp.LorentzianSusceptibility(frequency=cSi_frq1, gamma=cSi_gam1, sigma=cSi_sig1)]
cSi_susc += [mp.LorentzianSusceptibility(frequency=cSi_frq2, gamma=cSi_gam2, sigma=cSi_sig2)]
cSi_susc += [mp.LorentzianSusceptibility(frequency=cSi_frq3, gamma=cSi_gam3, sigma=cSi_sig3)]


med2 = mp.Medium(epsilon=1.0, E_susceptibilities=cSi_susc, valid_freq_range=cSi_range, \
    D_conductivity=Dcon/2/np.pi/2/np.pi)
# perm2 = np.array([med2.epsilon(1/w)[0,0] for w in wavem])
# perm2 = np.array([mp.Medium(epsilon=ceps.real, \
    # D_conductivity=2*np.pi/w * ceps.imag / ceps.real).epsilon(1/w)[0,0] for w in wavem])

nn = 3.859
kk = 0.015
ceps2 = (nn**2 - kk**2) + 1j*(2*nn*kk)


#Plot
fig, axes = plt.subplots(2, figsize=(6,9), sharex=True)
axes[0].plot(wave1, pexp1.real, '-',  label=ext1)
axes[0].plot(wave2, pexp2.real, '--', label=ext2)
axes[0].plot(wavem, perm1.real, '-.', label='Meep 1')
axes[0].plot(wavem, perm2.real, ':',  label='Meep 2')
axes[0].axhline(ceps2.real)
axes[0].axvline(0.641)
axes[0].set_ylabel(r'Re($\varepsilon)$')
axes[1].plot(wave1, pexp1.imag, '-')
axes[1].plot(wave2, pexp2.imag, '--')
axes[1].plot(wavem, perm1.imag, '-.')
axes[1].plot(wavem, perm2.imag, ':')
axes[1].axhline(ceps2.imag)
axes[1].axvline(0.641)
axes[1].set_ylabel(r'Im($\varepsilon)$')
axes[0].legend()


breakpoint()
