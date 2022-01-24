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
import meep as mp
import numpy as np
import matplotlib.pyplot as plt;plt.ion()
import h5py

class Analyzer(object):

    def __init__(self, params, prop=None, build_geo=False):
        self.util = semp.utils.Utilities()

        #Initialize
        self.initialize(params, prop=prop)

        #Get meta data
        if build_geo:
            self.build_geo()
        else:
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

        #Return the first index below wafer
        return np.where(self.xx >= obs_x)[0][0]

############################################
############################################

############################################
####	Process Data ####
############################################

    def get_data(self, comp, wave=None, ind=None, is_bbek=False):

        #X field has different vacuum data
        if comp[-1] == 'x':
            vac_comp = comp[0] + 'z'
        else:
            vac_comp = comp

        #Load data
        fld = self.load_field(comp, wave=wave, ind=ind)

        #Load vacuum data
        vac = self.load_field(vac_comp, wave=wave, ind=ind, is_vac=True)

        #Normalize by vacuum field
        fld /= vac

        #Subtract 1 if braunbek
        if is_bbek and comp[-1] != 'x':
            if self.prop.msim.geo.is_edge:
                fld -= np.heaviside(self.yy, 0)
            else:
                fld -= np.heaviside(self.yy, 0) * \
                    np.heaviside(self.prop.msim.gap_width - self.yy,0)

        return fld

    def get_fardata(self, comp, wave=None, is_bbek=False):
        #Load data
        fld = self.load_farfield(comp, wave=wave)
        vac = self.load_farfield(comp, wave=wave, is_vac=True)

        #Normalize by vacuum field
        if not np.allclose(np.abs(vac),0):
            fld /= vac

        #Subtract 1 if braunbek
        if is_bbek:
            if self.prop.msim.geo.is_edge:
                fld -= np.heaviside(self.far_yy, 0)
            else:
                fld -= np.heaviside(self.far_yy, 0) * \
                    np.heaviside(self.prop.msim.gap_width - self.far_yy,0)

        return fld

############################################
############################################

############################################
####	Load Data ####
############################################

    def load_metadata(self):

        #Load coordinates
        for is_vac in [True, False]:

            #Prefix
            pre = ['','vac-'][int(is_vac)]

            with h5py.File(f"{self.data_dir}/{pre}coords_waves.h5", 'r') as f:
                for c in ['xx','yy','zz']:
                    setattr(self, f"{['','vac_'][int(is_vac)]}{c}", f[c][()])
                self.waves = f['waves'][()]

        #Shift y to put zero at edge
        self.yy += self.prop.msim.geo.edge_y

        #Build farfield coordinates
        self.far_xx = self.prop.msim.farfield_z*self.util.m2mu
        self.far_yy = np.linspace(-0.5, 0.5, self.prop.msim.farfield_npts) * \
            self.prop.msim.farfield_width*self.util.m2mu + self.prop.msim.geo.edge_y

    ############################################

    def get_pol_wind(self, comp, wave, is_vac):

        #lowercase
        comp = comp.lower()

        #Vacuuum extension
        vac_ext = ['', 'vac-'][int(is_vac)]

        #polariztion
        if comp in ['ez', 'hx', 'hy']:
            pol = 's'
        else:
            pol = 'p'

        #Wave index
        wind = self.get_wind(wave)

        return pol, wind, comp, vac_ext

    def get_wind(self, wave=None):

        #wavelength index
        if wave is None:
            wind = 0
        else:
            wind = np.argmin(np.abs(wave - self.waves))
            #Check its close
            if np.abs(self.waves[wind] - wave)/wave > 0.01:
                print('\nWavelength is not close!\n')
                breakpoint()

        return wind

    ############################################

    def load_field(self, comp, wave=None, is_vac=False, ind=None):

        #Fix ind
        if ind is None:
            ind = slice(None)

        #Get polarization and wave index
        pol, wind, comp, vac_ext = self.get_pol_wind(comp, wave, is_vac)

        #Filename
        fname = f'{self.data_dir}/{vac_ext}fields_{pol}.h5'

        #Load data
        with h5py.File(fname, 'r') as f:
            data = f[f'{comp}_{wind}.r'][()] + 1j*f[f'{comp}_{wind}.i'][()]

        #Extract index slice
        data = data[ind]

        #Add shape to vacuum to divide by fld
        if is_vac and len(data.shape) != 0:
            data = data[:,None]

        return data

    ############################################

    def load_farfield(self, comp, wave=None, is_vac=False):

        #Get polarization and wave index
        pol, wind, comp, vac_ext = self.get_pol_wind(comp, wave, is_vac)

        #Filename
        fname = f'{self.data_dir}/{vac_ext}farfields_{pol}.h5'

        #Load data
        with h5py.File(fname, 'r') as f:
            data = f[f'{comp.capitalize()}'][:,wind]

        return data

############################################
############################################

############################################
####	Build Geometry ####
############################################

    def build_geo(self):

        #Build sim
        sim = self.prop.msim.build_sim()
        sim.init_sim()

        #Get dielectric
        self.eps = np.abs(sim.get_epsilon(self.prop.msim.fcen0))

        #Get coordinates
        self.xx, self.yy, self.zz, w = sim.get_array_metadata()

############################################
############################################

############################################
####	Collect Braunbek ####
############################################

    def collect_braunbek(self, wave=None):

        #Get xind at bottom of wafer (plus obs_distance)
        xind = self.get_xind()

        #Data to collect [[s-fld, s-drv], [p-fld, p-drv]]
        data_names = [['ez', 'hy'], ['hz', 'ey']]

        #Loop through polarizations
        data = []
        for data_list in data_names:

            #Get data
            fld = self.get_data(data_list[0], wave=wave, ind=xind, is_bbek=True)
            drv = self.get_data(data_list[1], wave=wave, ind=xind, is_bbek=True)

            avg = (fld + drv) / 2

            #Append
            data.append(avg)

        #Return s, p, yy (called x in diffraq)
        return data[0], data[1], self.yy

############################################
############################################

############################################
####	Plot Analyses ####
############################################

    def show_image(self, comp, wave=None, is_phase=False, is_bbek=False, vmax=None):

        #Load data
        data = self.get_data(comp, wave=wave, is_bbek=is_bbek)

        #Convert
        if is_phase:
            data = np.angle(data)
        else:
            data = np.abs(data)

        #Plot
        return data, self.plotter.plot_image(data, is_phase=is_phase, title=comp, vmax=vmax)

    def show_slice(self, comp, wave=None, is_phase=False, is_bbek=False):
        #Get index
        xind = self.get_xind()

        #Load data
        data = self.get_data(comp, wave=wave, ind=xind, is_bbek=is_bbek)

        #Convert
        if is_phase:
            data = np.angle(data)
        else:
            data = np.abs(data)

        #Plot
        return data, self.plotter.plot_slice(data, is_phase=is_phase)

    def show_epsilon(self, with_lines=True):

        #Check if epsilon already exists
        if not hasattr(self, 'eps'):
            self.build_geo()

        #Plot
        return self.plotter.plot_epsilon(self.eps, with_lines=with_lines)

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
        names = ['xx', 'yy', 'zz', 'eps']
        for nn in names:
            if hasattr(self, nn):
                delattr(self, nn)

############################################
############################################
