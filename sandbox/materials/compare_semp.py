import numpy as np
import semp
import matplotlib.pyplot as plt;plt.ion()

is_bbek = [False, True][1]
obs_x = 0

alz1 = semp.analysis.Analyzer({'session': 'materials/epsilon'})
xind = alz1.get_xind(obs_x=alz1.prop.msim.wafer_thick/2 + obs_x)
ez1 = alz1.get_data('ez', is_bbek=is_bbek)[xind]
hz1 = alz1.get_data('hz', is_bbek=is_bbek)[xind]

alz2 = semp.analysis.Analyzer({'session': 'materials/skin_meep_epsilon'})
ez2 = alz2.get_data('ez', is_bbek=is_bbek)[xind]
hz2 = alz2.get_data('hz', is_bbek=is_bbek)[xind]


fig, axes = plt.subplots(2, 2, figsize=(9,9))
axes[0,0].plot(alz1.yy, abs(ez1), '-',  label='Meep')
axes[0,0].plot(alz2.yy, abs(ez2), '--', label='Epsilon')

axes[1,0].plot(alz1.yy, abs(hz1), '-',  label='Meep')
axes[1,0].plot(alz2.yy, abs(hz2), '--', label='Epsilon')


axes[0,1].plot(alz1.yy, np.angle(ez1), '-',  label='Meep')
axes[0,1].plot(alz2.yy, np.angle(ez2), '--', label='Epsilon')

axes[1,1].plot(alz1.yy, np.angle(hz1), '-',  label='Meep')
axes[1,1].plot(alz2.yy, np.angle(hz2), '--', label='Epsilon')

axes[0,0].legend()

breakpoint()
