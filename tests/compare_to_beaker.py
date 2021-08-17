import numpy as np
import semp
import beaker
import matplotlib.pyplot as plt;plt.ion()

##SEMP
sparams = {
    'session':      'comp_M12P3/641nm',
}

salz = semp.analysis.Analyzer(sparams)
ssf, spf, sxx = salz.collect_braunbek()

##BEAKER
bparams = {
    'sim_type':         'braunbek',
    'polarization':     's',
    'session_name':     'Milestone_2/M12P3',
    'load_ext':         '641_ta5_sd225',
}

#Load data
balz = beaker.Analyzer(bparams)
balz.load_sim_instance()
bxx, alz_wav, bsf = balz.load_data_braunbek()
bparams['polarization'] = 'p'
balz = beaker.Analyzer(bparams)
balz.load_sim_instance()
bxx, alz_wav, bpf = balz.load_data_braunbek()

##Plot
fig, axes = plt.subplots(2, sharex=True, sharey=True, figsize=(6,9))
axes[0].plot(sxx, abs(ssf), label='s-SEMP')
axes[0].plot(bxx, abs(bsf), label='s-BEAKER')
axes[1].plot(sxx, abs(spf), label='p-SEMP')
axes[1].plot(bxx, abs(bpf), label='p-BEAKER')
axes[0].set_xlim([-5,10])
for i in range(2):
    axes[i].legend()
    axes[i].set_ylabel('Amplitude')

# fig2, axes2 = plt.subplots(2, sharex=True, sharey=True, figsize=(6,9))
# axes2[0].plot(sxx, np.angle(ssf), label='s-SEMP')
# axes2[0].plot(bxx, np.angle(bsf), label='s-BEAKER')
# axes2[1].plot(sxx, np.angle(spf), label='p-SEMP')
# axes2[1].plot(bxx, np.angle(bpf), label='p-BEAKER')
# for i in range(2):
#     axes2[i].legend()
#     axes2[i].set_ylabel('Phase')

#Integrate
sdx = sxx[1] - sxx[0]
bdx = bxx[1] - bxx[0]

# ssi = np.trapz(ssf[sxx > 0], dx=sdx)
# spi = np.trapz(spf[sxx > 0], dx=sdx)
# bsi = np.trapz(bsf[bxx > 0], dx=bdx)
# bpi = np.trapz(bpf[bxx > 0], dx=bdx)
ssi = np.trapz(ssf, dx=sdx)
spi = np.trapz(spf, dx=sdx)
bsi = np.trapz(bsf, dx=bdx)
bpi = np.trapz(bpf, dx=bdx)

print(abs(ssi), np.angle(ssi), abs(bsi), np.angle(bsi))
print(abs(spi), np.angle(spi), abs(bpi), np.angle(bpi))

breakpoint()
