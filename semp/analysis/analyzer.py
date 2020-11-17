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

class Analyzer(object):

    def __init__(self, params):
        self.util = semp.utils.Utilities()
        #Initialize
        self.initialize(params)

############################################
####	Initialization ####
############################################

    def initialize(self, params):
        #Set meep parameters
        self.util.set_default_params(self, params, semp.utils.def_params_ANLZ)
        self.params = params

        #Load directory
        self.data_dir = r"%s/%s"%(self.base_dir, self.session)

        #Load parameters (user-specified + default)
        prop_params = semp.utils.Logger.load_parameters(semp.utils.Logger, alz=self)

        #Create PROP instance
        self.prop = semp.Propagator(prop_params, is_analysis=True)

############################################
############################################

############################################
####	Shared Data functions ####
############################################

    def normalize_fields(self, pkg, is_difference=False):

        #Save metadata
        self.fld_comp = pkg['fld_comp']
        self.drv_comp = pkg['drv_comp']

        #Extract fields
        fld = pkg['waf_fld'].copy()
        drv = pkg['waf_drv'].copy()
        vac = pkg['vac_fld'][:, None]
        vdr = pkg['vac_drv'][:, None]

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

        #Normalize fields (absolute value, so -/+ for s/p is used in combining with main field)
        if fld.ndim == 3:
            vac = vac[:,None]
            vdr = vdr[:,None]
        fld /= np.abs(vac)
        drv /= np.abs(vdr)

        #Store relevant fields
        self.fld = fld
        self.drv = drv
        self.xx = pkg['waf_xx']
        self.yy = pkg['waf_yy']
        self.zz = pkg['waf_zz']

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
        self.normalize_fields(pkg, is_difference=False)

        #Grab slice (2 pixels wide) at observation point and average over
        dx = np.abs(self.xx - self.prop.msim.geo.obs_pt_x)
        xinds = np.isclose(dx, 0.5/self.prop.msim.resolution)
        self.fld = self.fld[xinds].mean(0)
        self.drv = self.drv[xinds].mean(0)
        self.xx = self.xx[xinds].mean(0)


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
