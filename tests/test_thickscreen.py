"""
test_thickscreen.py

Author: Anthony Harness
Affiliation: Princeton University
Created on: 11-19-2020
Package: SEMP

Description: Test SEMP package by simulating thick PEC screen and comparing to
    analytical solution.
License: Refer to $pkg_home_dir/LICENSE
"""

import numpy as np
import semp

class Test_Thick(object):

    is_debug = True

    ### HARDWIRED ###
    waves = [0.641, 0.725]
    pad = 4
    pml = 4.
    seam_dark = 4.
    seam_lite = 0.
    resolution = 40
    atol = 0.02
    wafer_thick = 2
    gap_width = 8
    do_plot = is_debug and semp.mpi_size == 1

############################################
####	Tests ####
############################################

    def test_all(self):

        #Run simulations (s & p)
        self.run_sim()

        #Loop through observation distances
        for obs_x in [0., -0.5]:

            #Get analyzer
            alz_params = {
                'base_dir':         f'{semp.tmp_dir}/tests',
                'session':          'thick_screen_noabs_1e-5',
                'obs_distance':     obs_x,
            }
            alz = semp.analysis.Analyzer(alz_params)

            #Shift yy to center
            alz.yy -= alz.prop.msim.geo.edge_y

            #Loop through waves
            for wave in self.waves[:]:

                #Loop through and check with and without braunbek
                for is_bbek in [False, True][:1]:

                    #Check analyzer
                    self.check_analyzer(wave, alz, is_bbek)

    ############################################

    def check_analyzer(self, wave, alz, is_bbek):

        #Get xindex
        xind = alz.get_xind()
        breakpoint()

        #Load normalized data
        data_list = ['ez','hz','ey','hy']
        sdata = [alz.get_data(dn, wave=wave, ind=xind, is_bbek=is_bbek) for dn in data_list]

        #Compare to sommerfeld
        adata, ddata = self.compare_to_analytic(wave, sdata, xind, alz, is_bbek)

        #Plot?
        if self.do_plot:
            fig, axes, fig2, axes2 = self.show_plot(sdata, adata, ddata, alz.yy)
            breakpoint()

    ############################################

    def compare_to_analytic(self, wave, sdata, xind, alz, is_bbek):

        #Load analytic
        thk = semp.analysis.Thick_Solution({'wave':wave,
            'wafer_thick':self.wafer_thick, 'gap_width':self.gap_width})

        #Get anaytic solution
        axx = alz.xx[xind]
        aex, aey, aez, ahx, ahy, ahz = thk.get_thick_solution(axx, alz.yy, \
            is_bbek=is_bbek)

        #Unpackage
        sez, shz, sey, shy = sdata

        #Get total diff
        dez = np.abs(sez - aez).mean()
        dhz = np.abs(shz - ahz).mean()
        dey = np.abs(sey - aey).mean()
        dhy = np.abs(shy - ahy).mean()

        #Assert all fields are true
        if not self.is_debug:
            assert(dez < self.atol and dhz < self.atol and dey < self.atol and dhy < self.atol)
        else:
            print(dez < self.atol and dhz < self.atol and dey < self.atol and dhy < self.atol)
            print(dez, dhz, dey, dhy)

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
            'sim_geometry':     'gap',
            'wafer_material':   'metal',
            'wafer_thick':      self.wafer_thick,
            'gap_width':        self.gap_width,
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
            'session':          'thick_screen',
        }

        #Run simulation
        prop = semp.Propagator(MEEP_params, PROP_params)
        if not self.is_debug:
            prop.run_sim()

        ##Debug Only
        if self.is_debug:
            if semp.mpi_size > 1:
                prop.run_sim()
                import sys;sys.exit(0)

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

        xlim = self.gap_width/2 + 1

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
        axes[0,0].set_xlim(-xlim, xlim)
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
        axes2[0,0].set_xlim(-xlim, xlim)
        for i in range(2):
            axes2[0,i].axhline(0.5, color='k', linestyle=':')
            axes2[0,i].legend()
            for j in range(2):
                axes2[j,i].axvline(0., color='k', linestyle=':')

        breakpoint()
        return fig, axes, fig2, axes2

############################################
############################################

if __name__ == '__main__':

    test = Test_Thick()
    test.test_all()
