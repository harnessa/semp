"""
analyzer.py

Author: Anthony Harness
Affiliation: Princeton University
Created on: 11-16-2020
Package: SEMP

Description: Class to load and analyze data from SEMP run
License: Refer to $pkg_home_dir/LICENSE
"""

import semp
import numpy as np
import matplotlib.pyplot as plt;plt.ion()
import h5py

class Analyzer(object):

    def __init__(self, params, prop=None):
        self.util = semp.utils.Utilities()

        #Initialize
        self.initialize(params, prop=prop)

        #Get meta data
        self.load_metadata()

############################################
####	Initialization ####
############################################

    def initialize(self, params, prop=None):
        #Set meep parameters
        self.util.set_default_params(self, params, semp.utils.def_params_ANLZ)
        self.params = params

        #Load directory
        self.data_dir = r"%s/%s"%(self.base_dir, self.session)

        #Store prop
        if prop is not None:
            self.prop = prop
        else:
            #Load parameters (user-specified + default)
            prop_params = semp.utils.Logger.load_parameters(semp.utils.Logger, alz=self)

            #Create PROP instance
            self.prop = semp.Propagator(prop_params, is_analysis=True)

        #Shift observation point to align with wafer bottom
        self.obs_distance += self.prop.msim.wafer_thick/2

        #Store data dir
        self.data_dir = self.prop.logger.data_dir

        #Save directory
        if self.do_save:
            if self.save_dir is None:
                self.save_dir = self.prop.logger.data_dir

        #Initialize plotter
        self.plotter = semp.analysis.Plotter(self)

############################################
############################################

############################################
####	Shared Functions ####
############################################

    def get_xind(self, obs_x=None):

        #Use default observation x if none supplied
        if obs_x is None:
            obs_x = self.obs_distance

        return np.argmin(np.abs(self.xx - 0.5/self.prop.msim.resolution - obs_x))

############################################
############################################

############################################
####	Process Data ####
############################################

    def get_data(self, comp, ind=None, is_bbek=False):
        #Load data
        fld = self.load_field(comp, ind=ind)
        vac = self.load_field(comp, ind=ind, is_vac=True)

        #Normalize by vacuum field
        fld /= vac

        #Turn subtract 1 if braunbek
        if is_bbek:
            fld -= np.heaviside(self.yy, 1)

        return fld

############################################
############################################

############################################
####	Load Data ####
############################################

    def load_metadata(self):

        #Get time extension for filenames
        time_ext = float(np.load(f'{self.data_dir}/time_ext.npy'))
        self.time_ext = f"-{time_ext:09.2f}"

        #Load coordinates
        for is_vac in [True, False]:
            with h5py.File(f"{self.data_dir}/{['','vac-'][int(is_vac)]}coords.h5", 'r') as f:
                for c in ['xx','yy','zz']:
                    setattr(self, f"{['','vac_'][int(is_vac)]}{c}", f[c][()])

        #Store pml size (with pad for y)
        self.pnum_x = int(self.prop.msim.geo.pmlx * self.prop.msim.resolution)
        self.pnum_y = int(self.prop.msim.geo.padpmly * self.prop.msim.resolution)
        self.pnum_z = int(self.prop.msim.geo.pmlz * self.prop.msim.resolution)

        #Trim pml
        self.xx = self.xx[self.pnum_x:self.xx.size-self.pnum_x]
        self.yy = self.yy[self.pnum_y:self.yy.size-self.pnum_y]
        self.zz = self.zz[self.pnum_z:self.zz.size-self.pnum_z]

        #Shift y to put zero at edge
        if self.prop.msim.sim_geometry == 'edge':
            self.yy += self.prop.msim.geo.edge_y

    def load_field(self, comp, is_vac=False, ind=None):

        #Fix ind
        if ind is None:
            ind = slice(None)

        #lowercase
        comp = comp.lower()

        #Vacuuum extension
        vac_ext = ['', 'vac-'][int(is_vac)]

        #Filename
        fname = self.data_dir + '/' + vac_ext + comp + self.time_ext + ".h5"

        #Load data (without pml) - ugly index is due to slow fancy indexing in h5py
        with h5py.File(fname, 'r') as f:
            if not is_vac:
                sx, sy = f[f'{comp}.r'].shape
                data = f[f'{comp}.r'][self.pnum_x:sx-self.pnum_x, self.pnum_y:sy-self.pnum_y] + \
                    1j*f[f'{comp}.i'][self.pnum_x:sx-self.pnum_x, self.pnum_y:sy-self.pnum_y]
            else:
                sx = f[f'{comp}.r'].shape[0]
                data = f[f'{comp}.r'][self.pnum_x:sx-self.pnum_x] + \
                    1j*f[f'{comp}.i'][self.pnum_x:sx-self.pnum_x]

        #Extract index slice
        data = data[ind]

        #Add shape to vacuum to divide by fld
        if is_vac and len(data.shape) != 0:
            data = data[:,None]

        return data

############################################
############################################

############################################
####	Plot Analyses ####
############################################

    def show_image(self, comp, is_phase=False, is_bbek=False):

        #Load data
        data = self.get_data(comp, is_bbek=is_bbek)

        #Convert
        if is_phase:
            data = np.angle(data)
        else:
            data = np.abs(data)

        #Plot
        return self.plotter.plot_image(data, is_phase=is_phase)

    def show_slice(self, comp, is_phase=False, is_bbek=False):
        #Get index
        xind = self.get_xind()

        #Load data
        data = self.get_data(comp, ind=xind, is_bbek=is_bbek)

        #Convert
        if is_phase:
            data = np.angle(data)
        else:
            data = np.abs(data)

        #Plot
        return self.plotter.plot_slice(data, is_phase=is_phase)

############################################
############################################

############################################
####	Movie Analysis ####
############################################

    def analyze_movie(self):
        pass

############################################
############################################

############################################
####	Misc ####
############################################

    def clean_up(self):
        names = ['xx', 'yy', 'zz']
        for nn in names:
            if hasattr(self, nn):
                delattr(self, nn)

############################################
############################################
