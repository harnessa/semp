import semp
import numpy as np
import matplotlib.pyplot as plt;plt.ion()

wave = [0.641, 0.725][1]
obsx = 0

#Load Semp params
alz_params = {
    'base_dir':         f'{semp.tmp_dir}/tests',
    'session':          'thick_screen_all',
    'obs_distance':     obsx,
}
alz = semp.analysis.Analyzer(alz_params)

gap_width = alz.prop.msim.gap_width
wafer_thick = alz.prop.msim.wafer_thick

#Load metadata
xind = alz.get_xind()
xx = alz.xx[xind]
yy = alz.yy.copy()
yy -= alz.prop.msim.geo.edge_y

#Trim to gap
ginds = np.abs(yy) <= gap_width/2.
yy = yy[ginds]
yy_plt = yy.copy()
yy += gap_width/2

#Load data
if [False, True][1]:
    sez = alz.get_data('ez', wave=wave, ind=xind)[ginds]
    shy = alz.get_data('hy', wave=wave, ind=xind)[ginds]
    #Need to normalize by hx by hy vacuum
    shx, _ = alz.load_field('hx', wave=wave, ind=xind, is_vac=False)
    vhy, _ = alz.load_field('hy', wave=wave, ind=xind, is_vac=True)
    shx = shx[ginds] / vhy

else:
    sez = None

# bez = alz.get_data('ez', wave=wave)
# plt.imshow(abs(bez))
# breakpoint()

########################################
########################################


#Shared calcs
zz = np.abs(xx) + wafer_thick/2.
kk = 2.*np.pi/wave
H0 = 1
E0 = 1
ky = lambda n: np.pi*n/gap_width

beta = lambda n: np.sqrt(kk**2. - ky(n)**2.)
cos_n = lambda n: np.cos(np.pi*n) - 1.
exp_n = lambda n: np.exp(1j*beta(n)*zz)
sqrt_mu_eps = 1.
om_mu = kk * sqrt_mu_eps
om_ep = kk / sqrt_mu_eps

################

#TE (s-pol)
Bn = lambda n: -1j*2./(gap_width*beta(n)) * H0 * cos_n(n)
Ez_func1 = lambda n:  1j*om_mu / ky(n) * Bn(n) * np.sin(ky(n)*yy) * exp_n(n)
Hy_func1 = lambda n: -1j*beta(n) / ky(n) * Bn(n) * np.sin(ky(n)*yy) * exp_n(n)
Hx_func1 = lambda n: Bn(n) * np.cos(ky(n)*yy) * exp_n(n)

#Works
# Ez_func2 = lambda n:  4*E0/(n*np.pi) * np.sin(ky(n)*yy) * exp_n(n) * (kk / beta(n))
# Hy_func2 = lambda n: -4*E0/(n*np.pi) * np.sin(ky(n)*yy) * exp_n(n)
# Hx_func2 = lambda n: -4j*E0/(gap_width*om_mu) * np.cos(ky(n)*yy) * exp_n(n) * (kk / beta(n))

#Works best with correct signs too
E0 = lambda n: om_mu / beta(n) * H0
Ez_const = lambda n: 4*E0(n)/(n*np.pi)
Ez_func2 = lambda n: Ez_const(n) * np.sin(ky(n)*yy) * exp_n(n)
Hy_func2 = lambda n: beta(n)/om_mu * Ez_func2(n)
Hx_func2 = lambda n: 1j*ky(n)/om_mu * Ez_const(n) * np.cos(ky(n)*yy) * exp_n(n)

##hx doesnt work
# Ez_func2 = lambda n: 4/(n*np.pi) * (om_mu/beta(n)) * np.sin(ky(n)*yy) * np.exp(1j*beta(n)*zz)
# Hy_func2 = lambda n: 4/(n*np.pi) * np.sin(ky(n)*yy) * np.exp(1j*beta(n)*zz)
# Hx_func2 = lambda n: -4j/(n*np.pi) * (ky(n)**2/beta(n)) * np.cos(ky(n)*yy) * np.exp(1j*beta(n)*zz)

# works
# E0 = lambda n: (kk / beta(n))
# Ez_func2 = lambda n:  4*E0(n)/(gap_width*ky(n)) * np.sin(ky(n)*yy) * np.exp(1j*beta(n)*zz)
# Hy_func2 = lambda n: -4*E0(n)*beta(n)/(gap_width*ky(n)*om_mu) * np.sin(ky(n)*yy) * np.exp(1j*beta(n)*zz)
# Hx_func2 = lambda n: 4j*E0(n)/(gap_width*om_mu) * np.cos(ky(n)*yy) * np.exp(1j*beta(n)*zz)

#Get mode numbers
Nn = int(np.ceil(2.*gap_width/wave))
nns = np.arange(1,Nn, 2)            #even n's are zero b/c of cos_n(n)
nns = nns.reshape(len(nns), 1)

#Calculate fields
ez1 = Ez_func1(nns).sum(0)
hx1 = Hx_func1(nns).sum(0)
hy1 = Hy_func1(nns).sum(0)

ez2 = Ez_func2(nns).sum(0)
hx2 = Hx_func2(nns).sum(0)
hy2 = Hy_func2(nns).sum(0)

#Normalize by plane wave
pw1 = np.exp(1j*kk*zz)
ez1 /= pw1
hy1 /= pw1
hx1 /= pw1
pw2 = np.exp(1j*kk*zz)
ez2 /= pw2
hy2 /= pw2
hx2 /= pw2


#FIXME: why are some of these flipped?
ez1 *= -1
# hy2 *= -1
# hx2 *= -1


#######################################
#######################################

fig, axes = plt.subplots(3, 2, figsize=(8,11))

axes[0,0].plot(abs(ez1))
axes[0,0].plot(abs(ez2), '--')
axes[0,0].set_title('Ez Amp')

axes[1,0].plot(abs(hy1))
axes[1,0].plot(abs(hy2), '--')
axes[1,0].set_title('Hy Amp')

axes[2,0].plot(abs(hx1))
axes[2,0].plot(abs(hx2), '--')
axes[2,0].set_title('Hx Amp')

axes[0,1].plot(np.angle(ez1))
axes[0,1].plot(np.angle(ez2), '--')
axes[0,1].set_title('Ez Phase')

axes[1,1].plot(np.angle(hy1))
axes[1,1].plot(np.angle(hy2), '--')
axes[1,1].set_title('Hy Phase')

axes[2,1].plot(np.angle(hx1))
axes[2,1].plot(np.angle(hx2), '--')
axes[2,1].set_title('Hx Phase')

if sez is not None:
    axes[0,0].plot(abs(sez), ':')
    axes[1,0].plot(abs(shy), ':')
    axes[2,0].plot(abs(shx), ':')
    axes[0,1].plot(np.angle(sez), ':')
    axes[1,1].plot(np.angle(shy), ':')
    axes[2,1].plot(np.angle(shx), ':')

breakpoint()
