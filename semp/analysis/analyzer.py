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

        #Store data dir
        self.data_dir = self.prop.logger.data_dir

        #Save directory
        if self.do_save:
            if self.save_dir is None:
                self.save_dir = self.prop.logger.data_dir

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

        #Shift if edge
        if self.prop.msim.sim_geometry == 'edge':
            self.yy += self.prop.msim.geo.edge_y

############################################
############################################

############################################
####	Process Data ####
############################################

    def get_data_slices(self, data_list, obs_x=None, is_bbek=False):

        #Use default observation x if none supplied
        if obs_x is None:
            obs_x = self.obs_distance

        #Find index corresponding to observation distance (account for Yee lattice)
        xind = np.argmin(np.abs(self.xx - 0.5/self.prop.msim.resolution - obs_x))

        #Load data
        data = []
        for dn in data_list:
            #Load data and normalize by vacuum field
            fld = self.load_field(dn, ind=xind) / \
                self.load_field(dn, ind=xind, is_vac=True)

            #Turn subtract 1 if braunbek
            if is_bbek:
                fld -= np.heaviside(self.yy, 1)

            #Append
            data.append(fld)

        return data, xind

############################################
############################################

############################################
####	Load Data ####
############################################

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

        return data

############################################
############################################

############################################
####	Normalize Field ####
############################################

    def normalize_fields(self, pkg, is_difference=False):

        #Save metadata
        self.fld_comp = pkg['fld_comp']
        self.drv_comp = pkg['drv_comp']

        #Extract fields
        fld = pkg['waf_fld'].copy()
        drv = pkg['waf_drv'].copy()
        vac = pkg['vac_fld']
        vdr = pkg['vac_drv']

        #Create difference fields
        if is_difference:

            #Build geometric shadow
            shadow = pkg['waf_eps'].copy()
            shadow[shadow != 1.] = 0.

            #Compress along x-axis (meep-x) to see where there isn't air all the way to the source
            shadow = shadow.sum(0)

            #Get where in geometric shadow (each pixel = 1, vertically)
            in_light = shadow == pkg['waf_xx'].shape[0]

            #Subtract out vacuum fields in air
            fld[..., in_light] -= vac
            drv[..., in_light] -= vdr

        #Normalize fields
        if fld.ndim == 3:
            vac = vac[:,None]
            vdr = vdr[:,None]

        fld /= vac#[:,None]
        drv /= vdr#[:,None]

        #Store relevant fields
        self.fld = fld
        self.drv = drv
        self.xx = pkg['waf_xx']
        self.yy = pkg['waf_yy']
        self.zz = pkg['waf_zz']
        self.eps = pkg['waf_eps']

        #Cleanup
        del pkg, vac, vdr

############################################
############################################

############################################
####	Difference Field Analysis ####
############################################

    def export_difference_field(self):
        #Load data
        pkg = self.prop.logger.load_propagation_data()

        #Normalize data
        self.normalize_fields(pkg, is_difference=True)

        #Get average field for Braunbek
        fld = (self.fld + self.drv)/2

        #Save data
        if self.do_save:
            #TODO: This could be redone, but naming is kept to be consistent with BEAKER
            with h5py.File(f'{self.save_dir}/edge{self.save_ext}.h5', 'w') as f:
                f.create_dataset('xx', data=self.yy)    #Edge distance
                f.create_dataset('fld', data=fld)
                f.create_dataset('wav', data=self.prop.msim.wave)
                f.create_dataset('zz', data=self.zz)

        plt.plot(self.yy, np.abs(self.fld))
        plt.plot(self.yy, np.abs(self.drv))
        import pdb;pdb.set_trace()

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
        names = ['xx', 'yy', 'zz', 'fld', 'drv']
        for nn in names:
            if hasattr(self, nn):
                delattr(self, nn)

############################################
############################################
