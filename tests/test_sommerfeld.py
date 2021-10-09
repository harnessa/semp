"""
test_sommerfeld.py

Author: Anthony Harness
Affiliation: Princeton University
Created on: 11-19-2020
Package: SEMP

Description: Test SEMP package by simulating thin PEC screen and comparing to
    Sommerfeld's analytical solution.
License: Refer to $pkg_home_dir/LICENSE
"""

import numpy as np
import semp

class Test_Sommerfeld(object):

    ### HARDWIRED ###
    waves = [0.641, 0.725]
    pad = 4.
    pml = 4.
    seam_dark = 5.
    seam_lite = 10.
    resolution = 30
    atol = 0.02
    do_plot = False and semp.mpi_size == 1

############################################
####	Tests ####
############################################

    def test_all(self):

        #Run simulations (s & p)
        self.run_sim()

        #Loop through observation distances
        for obs_x in [0, 3]:

            #Get analyzer
            alz_params = {
                'base_dir':         f'{semp.tmp_dir}/tests',
                'session':          'sommer',
                'obs_distance':     obs_x,
            }
            alz = semp.analysis.Analyzer(alz_params)

            #Loop through waves
            for wave in self.waves:

                #Loop through and check with and without braunbek
                for is_bbek in [False, True]:

                    #Check simulation
                    self.check_simulation(wave, alz, is_bbek)

                    #Check analyzer
                    self.check_analyzer(wave, alz, is_bbek)

    def check_simulation(self, wave, alz, is_bbek):

        #Find index corresponding to observation distance (account for Yee lattice)
        xind = np.argmin(np.abs(alz.xx - 0.5/alz.prop.msim.resolution - alz.obs_distance))

        #Load and normalize data
        sez = alz.load_field('ez', wave=wave, ind=xind) / \
            alz.load_field('ez', wave=wave, ind=xind, is_vac=True)
        shz = alz.load_field('hz', wave=wave, ind=xind) / \
            alz.load_field('hz', wave=wave, ind=xind, is_vac=True)
        sey = alz.load_field('ey', wave=wave, ind=xind) / \
            alz.load_field('ey', wave=wave, ind=xind, is_vac=True)
        shy = alz.load_field('hy', wave=wave, ind=xind) / \
            alz.load_field('hy', wave=wave, ind=xind, is_vac=True)

        #Subtract Braunbek field
        if is_bbek:
            sez -= np.heaviside(alz.yy, 1)
            shz -= np.heaviside(alz.yy, 1)
            sey -= np.heaviside(alz.yy, 1)
            shy -= np.heaviside(alz.yy, 1)

        #Package
        sdata = [sez, shz, sey, shy]

        #Compare to sommerfeld
        adata, ddata = self.compare_to_sommerfeld(wave, sdata, xind, alz, is_bbek)

        #Plot?
        if self.do_plot:
            fig, axes, fig2, axes2 = self.show_plot(sdata, adata, ddata, alz.yy)
            breakpoint()

    ############################################

    def check_analyzer(self, wave, alz, is_bbek):

        #Get xindex
        xind = alz.get_xind()

        #Load normalized data
        data_list = ['ez','hz','ey','hy']
        sdata = [alz.get_data(dn, wave=wave, ind=xind, is_bbek=is_bbek) for dn in data_list]

        #Compare to sommerfeld
        adata, ddata = self.compare_to_sommerfeld(wave, sdata, xind, alz, is_bbek)

        #Plot?
        if self.do_plot:
            fig, axes, fig2, axes2 = self.show_plot(sdata, adata, ddata, alz.yy)
            breakpoint()

    ############################################

    def compare_to_sommerfeld(self, wave, sdata, xind, alz, is_bbek):

        #Load sommerfeld
        som = semp.analysis.Sommerfeld({'wave':wave})

        #Get sommerfeld solution
        axx = alz.xx[xind]

        aex, aey, aez, ahx, ahy, ahz = som.get_sommerfeld_solution(axx, alz.yy, \
            is_bbek=is_bbek)

        #Unpackage
        sez, shz, sey, shy = sdata

        #Get total diff
        dez = np.abs(sez - aez).mean()
        dhz = np.abs(shz - ahz).mean()
        dey = np.abs(sey - aey).mean()
        dhy = np.abs(shy - ahy).mean()

        #Assert all fields are true
        assert(dez < self.atol and dhz < self.atol and dey < self.atol and dhy < self.atol)

        #Package
        adata = [aez, ahz, aey, ahy]
        ddata = [dez, dhz, dey, dhy]

        return adata, ddata

############################################
############################################

############################################
####	Numerical Solution ####
############################################

    def run_sim(self):
        MEEP_params = {
            ### Lab Properties  ###
            'polars':           ['s', 'p'],
            'waves':            self.waves,

            ### Mask Properties ###
            'sim_geometry':     'edge',
            'is_sommerfeld':    True,
            'seam_dark':        self.seam_dark,
            'seam_lite':        self.seam_lite,

            ### Numerics ###
            'resolution':       self.resolution,
            'pml_all':          self.pml,
            'pad_all':          self.pad,
            'decay_dt':         30,
            'use_absorber':     True,
        }

        PROP_params = {
            'verbose':          True,
            'base_dir':         f'{semp.tmp_dir}/tests',
            'session':          'sommer',
        }

        #Run simulation
        prop = semp.Propagator(MEEP_params, PROP_params)
        prop.run_sim()

        # ##Debug Only
        # if semp.mpi_size > 1:
        #     prop.run_sim()
        #     import sys;sys.exit(0)

############################################
############################################

############################################
####	Plotting ####
############################################

    def show_plot(self, sdata, adata, ddata, yy):

        import matplotlib.pyplot as plt;plt.ion()

        #Get data out
        sez, shz, sey, shy = sdata
        aez, ahz, aey, ahy = adata
        dez, dhz, dey, dhy = ddata

        #Print
        print(f'\n{dez:.2e}, {dhz:.2e}, {dey:.2e}, {dhy:.2e}\n')
        print(dez < self.atol and dhz < self.atol and dey < self.atol and dhy < self.atol)

        #Plot
        fig, axes = plt.subplots(2, 2, figsize=(9,9), sharex=True)
        axes[0,0].plot(yy, np.abs(aez), '-',  label='Sommer Ez')
        axes[0,0].plot(yy, np.abs(sez), '--',  label='Semp Ez')
        axes[1,0].plot(yy, np.angle(aez), '-',  label='Sommer Ez')
        axes[1,0].plot(yy, np.angle(sez), '--',  label='Semp Ez')

        axes[0,1].plot(yy, np.abs(ahz), '-',  label='Sommer Hz')
        axes[0,1].plot(yy, np.abs(shz), '--',  label='Semp Hz')
        axes[1,1].plot(yy, np.angle(ahz), '-',  label='Sommer Hz')
        axes[1,1].plot(yy, np.angle(shz), '--',  label='Semp Hz')

        axes[0,0].set_ylabel('Amplitude')
        axes[1,0].set_ylabel('Phase')
        axes[0,0].set_xlim(-3, 3)
        for i in range(2):
            axes[0,i].axhline(0.5, color='k', linestyle=':')
            axes[0,i].legend()
            for j in range(2):
                axes[j,i].axvline(0., color='k', linestyle=':')

        fig2, axes2 = plt.subplots(2, 2, figsize=(9,9), sharex=True)
        axes2[0,0].plot(yy, np.abs(aey), '-',  label='Sommer Ey')
        axes2[0,0].plot(yy, np.abs(sey), '--',  label='Semp Ey')
        axes2[1,0].plot(yy, np.angle(aey), '-',  label='Sommer Ey')
        axes2[1,0].plot(yy, np.angle(sey), '--',  label='Semp Ey')

        axes2[0,1].plot(yy, np.abs(ahy), '-',  label='Sommer Hy')
        axes2[0,1].plot(yy, np.abs(shy), '--',  label='Semp Hy')
        axes2[1,1].plot(yy, np.angle(ahy), '-',  label='Sommer Hy')
        axes2[1,1].plot(yy, np.angle(shy), '--',  label='Semp Hy')

        axes2[0,0].set_ylabel('Amplitude')
        axes2[1,0].set_ylabel('Phase')
        axes2[0,0].set_xlim(-3, 3)
        for i in range(2):
            axes2[0,i].axhline(0.5, color='k', linestyle=':')
            axes2[0,i].legend()
            for j in range(2):
                axes2[j,i].axvline(0., color='k', linestyle=':')

        return fig, axes, fig2, axes2

############################################
############################################

if __name__ == '__main__':

    test = Test_Sommerfeld()
    test.test_all()
