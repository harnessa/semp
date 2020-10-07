"""
propagator.py

Author: Anthony Harness
Affiliation: Princeton University
Created on: 10-06-2020
Package: SEMP

Description: Main class to control the entire SEMP python package
License: Refer to $pkg_home_dir/LICENSE
"""

import semp
import numpy as np
import atexit

class Propagator(object):

    def __init__(self, meep_params, prop_params=None, is_analysis=False):
        self.util = semp.utils.util
        self.is_analysis = is_analysis
        self.combine_params(meep_params, prop_params)
        self.initialize()

############################################
####	Initialize  ####
############################################

    def combine_params(self, meep_params, prop_params):
        if prop_params is not None:
            params = {'MEEP_params':meep_params, 'PROP_params':prop_params}
        else:
            params = meep_params
        self.params = self.util.deepcopy(params)

    def initialize(self):
        #Set propagation parameters
        self.util.set_default_params(self, \
            self.params['PROP_params'], semp.utils.def_params['PROP_params'])

        #Load Meep simulation class
        self.meep_sim = semp.simulation.Meep_Sim(self, self.params['MEEP_params'])

        #Load logger class
        self.logger = semp.utils.Logger(self)

############################################
############################################

############################################
####	Main Scripts ####
############################################

    def run_sim(self):
        #Run Movie or Propagation simulation
        if self.is_movie:
            self.run_movie()
        else:
            self.run_propagation()

    def get_epsilon(self):
        sim = self.meep_sim.build_sim()
        sim.init_fields()
        return sim.get_epsilon(self.meep_sim.fcen)

############################################
############################################

############################################
####	Propagation Simulation ####
############################################

    def run_propagation(self):
        import pdb;pdb.set_trace()

############################################
############################################

############################################
####	Movie Simulation ####
############################################

    def run_movie(self):
        #Run sim in vacuum
        self.run_single_movie(is_vac=True)

        #Run sim with wafer
        self.run_single_movie(is_vac=False)

    def run_single_movie(self, is_vac=False):
        #Build simulation
        sim = self.meep_sim.build_sim(is_vac=is_vac)

        #Set output directory + prefix
        sim.use_output_directory(self.save_dir)
        sim.filename_prefix = ''

        #Get filename
        filename = "%smovie_%s%s"%(['','vac_'][is_vac], \
            self.fld_comp_name, self.util.get_ext(self.save_ext))

        #Output function to save data
        if self.polarization == 's':
            fld_func = mp.output_efield_z
        else:
            fld_func = mp.output_hfield_z

        out_funcs = [mp.to_appended(filename, mp.at_every(self.meep_sim.save_dt, fld_func))]

        #Run sim
        sim.run(*out_funcs, until=self.meep_sim.run_time)

        #Reset sim
        sim.reset_meep()

        import pdb;pdb.set_trace()

############################################
############################################
