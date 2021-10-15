import numpy as np
import matplotlib.pyplot as plt;plt.ion()
import semp

#Get analyzer
alz_params = {
    'base_dir':         f'{semp.tmp_dir}/tests',
    'session':          'sommer',
    'obs_distance':     0,
}
alz = semp.analysis.Analyzer(alz_params)

wave = 0.641
kk = 2*np.pi/wave

#Load and normalize data
wez = alz.load_field('ez', wave=wave)
vez = alz.load_field('ez', wave=wave, is_vac=True)
whz = alz.load_field('hz', wave=wave)
vhz = alz.load_field('hz', wave=wave, is_vac=True)
wey = alz.load_field('ey', wave=wave)
vey = alz.load_field('ey', wave=wave, is_vac=True)
why = alz.load_field('hy', wave=wave)
vhy = alz.load_field('hy', wave=wave, is_vac=True)

sez = np.sign(vez[:,0])
shz = np.sign(vhz[:,0])
sey = np.sign(vey[:,0])
shy = np.sign(vhy[:,0])

#Normalize
if [False, True][0]:

    wez /= vez
    whz /= vhz
    wey /= vey
    why /= vhy

    # wez -= 1
    # whz -= 1
    # wey -= 1
    # why -= 1

#Get derivatives
dhy = -np.gradient(vez[:,0], alz.xx) / (1j*kk)
dey =  np.gradient(vhz[:,0], alz.xx) / (1j*kk)

fig, axes = plt.subplots(2, 2, figsize=(8,8))
axes[0,0].plot(alz.xx, abs(vhy), '-',  label='Hy - Meep')
axes[0,0].plot(alz.xx, abs(dhy), '--', label='Hy - Derv')
axes[1,0].plot(alz.xx, abs(vey), '-',  label='Ey - Meep')
axes[1,0].plot(alz.xx, abs(dey), '--', label='Ey - Derv')

axes[0,1].plot(alz.xx, np.angle(vhy), '-',  label='Hy - Meep')
axes[0,1].plot(alz.xx, np.angle(dhy), '--', label='Hy - Derv')
axes[1,1].plot(alz.xx, np.angle(vey), '-',  label='Ey - Meep')
axes[1,1].plot(alz.xx, np.angle(dey), '--', label='Ey - Derv')

axes[0,0].legend()
axes[1,0].legend()


#Get derivatives
# dhy2 = (np.gradient(wez[:,400], alz.xx)*vez[:,0] - np.gradient(vez[:,0], alz.xx)*wez[:,400])/vez[:,0]**2
# dhy2 *= -1./(1j*kk)
dhy2 = -np.gradient(wez[:,400], alz.xx) / (1j*kk)
dey2 =  np.gradient(whz[:,400], alz.xx) / (1j*kk)

if [False, True][0]:
    wez /= vez
    whz /= vhz
    wey /= vey
    why /= vhy
    dhy2 /= vhy[:,0]
    dey2 /= vey[:,0]

fig, axes = plt.subplots(2, 2, figsize=(8,8))
axes[0,0].plot(alz.xx, abs(why)[:,400], '-',  label='Hy - Meep')
axes[0,0].plot(alz.xx, abs(dhy2), '--', label='Hy - Derv')
axes[1,0].plot(alz.xx, abs(wey)[:,400], '-',  label='Ey - Meep')
axes[1,0].plot(alz.xx, abs(dey2), '--', label='Ey - Derv')

axes[0,1].plot(alz.xx, np.angle(why)[:,400], '-',  label='Hy - Meep')
axes[0,1].plot(alz.xx, np.angle(dhy2), '--', label='Hy - Derv')
axes[1,1].plot(alz.xx, np.angle(wey)[:,400], '-',  label='Ey - Meep')
axes[1,1].plot(alz.xx, np.angle(dey2), '--', label='Ey - Derv')

axes[0,0].legend()
axes[1,0].legend()

fig, axes = plt.subplots(2, 2, figsize=(8,8))
axes[0,0].plot(alz.xx, abs(wez[:,400]))
axes[0,0].plot(alz.xx, abs(why[:,400]))
axes[0,0].plot(alz.xx, abs(wez[:,400] - why[:,400])/2, '--')
axes[1,0].plot(alz.xx, abs(whz[:,400]))
axes[1,0].plot(alz.xx, abs(wey[:,400]))
axes[1,0].plot(alz.xx, abs(whz[:,400] + wey[:,400])/2, '--')

axes[0,1].plot(alz.xx, np.angle(wez[:,400]))
axes[0,1].plot(alz.xx, np.angle(why[:,400]))
axes[0,1].plot(alz.xx, np.angle(wez[:,400] - why[:,400])/2, '--')
axes[1,1].plot(alz.xx, np.angle(whz[:,400]))
axes[1,1].plot(alz.xx, np.angle(wey[:,400]))
axes[1,1].plot(alz.xx, np.angle(whz[:,400] + wey[:,400])/2, '--')

breakpoint()
