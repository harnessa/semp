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
import meep as mp

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
        self.msim = semp.simulation.Meep_Sim(self, self.params['MEEP_params'])

        #Load logger class
        self.logger = semp.utils.Logger(self)

        #Verbosity
        if not self.verbose:
            mp.verbosity(0)

############################################
############################################

############################################
####	Main Scripts ####
############################################

    def run_sim(self):
        #Run Startups
        self.logger.start_up()

        #Run Movie or Propagation simulation
        if self.is_movie:
            self.run_movie()
        else:
            self.run_propagation()

        #Run Closeups
        self.logger.close_up()

    def get_epsilon(self):
        sim = self.msim.build_sim()
        sim.init_fields()
        return sim.get_epsilon(self.msim.fcen)

############################################
############################################

############################################
####	Propagation Simulation ####
############################################

    def run_propagation(self):

        #Run sim in vacuum
        vac_data = self.run_single_propagation(is_vac=True)

        #Run sim with wafer
        waf_data = self.run_single_propagation(is_vac=False)

        #Data package
        data_pkg = {'vac':vac_data, 'waf':waf_data}

        #Save data
        self.logger.save_propagation_data(data_pkg)

        #Store data
        if self.store_results:
            self.data_pkg = self.logger.convert_data_package(data_pkg)
        else:
            #Cleanup
            del vac_data, waf_data, data_pkg

    def run_single_propagation(self, is_vac=False):

        #Build simulation
        sim = self.msim.build_sim(is_vac=is_vac)

        #Run sim
        sim.run(until=self.msim.run_time)

        #Synchronize magnetic field
        sim.fields.synchronize_magnetic_fields()

        #Get ouput volume
        vol = self.msim.geo.get_output_volume(is_vac=is_vac)

        #Get field array slice
        fld = sim.get_array(vol=vol, component=self.msim.src_comp, cmplx=True)

        #Get derivative field array slice
        drv = sim.get_array(vol=vol, component=self.msim.drv_comp, cmplx=True)

        #Get observation points (don't return w)
        xx, yy, zz, ww = sim.get_array_metadata(vol=vol)

        #Adjust coordinates
        xx, yy, zz = self.msim.adjust_coordinates(xx, yy, zz)

        #Get material map (epsilon)
        eps = sim.get_array(vol=vol, component=mp.Dielectric)

        #Collapse over y/z in vacuum
        if is_vac:
            #Double check it is the same
            if fld.std(1).max() != 0:
                print('\nVacuum not uniform over Y!\n')
            fld = fld.mean(1)
            drv = drv.mean(1)
            yy = np.mean(yy)
            zz = np.mean(zz)

        #Reset sim
        sim.reset_meep()

        #Cleanup
        del sim, ww

        #Build return package
        ret_pkg = {'fld':fld, 'drv':drv, 'eps':eps, 'xx':xx, 'yy':yy, 'zz':zz}

        return ret_pkg

############################################
############################################

############################################
####	Movie Simulation ####
############################################

    def run_movie(self):
        #Load movie maker
        self.movie_maker = semp.analysis.Movie_Maker(self)

        #Run sim with wafer
        self.run_single_movie()

    def run_single_movie(self):
        #Is vacuum?
        is_vac = self.msim.sim_geometry == 'vacuum'

        #Build simulation
        sim = self.msim.build_sim(is_vac=is_vac)

        #Set output directory + prefix
        sim.use_output_directory(self.logger.data_dir)
        sim.filename_prefix = ''

        #Get filename
        vac_ext = ['','vac_'][is_vac]
        filename = f'{vac_ext}movie_{self.msim.src_comp_name}{self.logger.log_ext}'

        #Function to output epsilon
        eps_func = mp.at_beginning(mp.output_epsilon)

        #Function to save field (appended) to h5 file
        if self.msim.polarization == 's':
            fld_func = mp.output_efield_z
        else:
            fld_func = mp.output_hfield_z
        h5_func = mp.to_appended(filename, mp.at_every(self.msim.save_dt, fld_func))

        #Function to save png files
        png_opts = "-Zc dkbluered -a yarg -C $EPS"
        png_func = mp.at_every(self.msim.save_dt, \
            mp.output_png(self.msim.src_comp, png_opts, rm_h5=True))

        #Output function to save data
        out_funcs = [eps_func, h5_func, png_func]

        #Run sim
        sim.run(*out_funcs, until=self.msim.run_time)

        #Reset sim
        sim.reset_meep()

        #Create movie
        self.movie_maker.make_movie(filename)

############################################
############################################
