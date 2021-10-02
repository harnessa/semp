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
import h5py

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
            self.run_to_end()

        #Run Closeups
        self.logger.close_up()

    def get_epsilon(self, sim=None):
        if sim is None:
            sim = self.msim.build_sim()
            sim.init_sim()
        return sim.get_epsilon(self.msim.fcen)

############################################
############################################

############################################
####	Run to End of Simulation ####
############################################

    def run_to_end(self):

        #Flag to get metadata
        get_meta = True

        #Loop over polarizations
        for pol in self.msim.polars:

            #Loop over vacuum
            for is_vac in [True, False]:

                #Run sim
                self.run_single_to_end(pol, is_vac, get_meta)

            #Turn off flag
            get_meta = False

    def run_single_to_end(self, pol, is_vac, get_meta):

        #Build to simulation
        sim = self.msim.build_sim(pol=pol, is_vac=is_vac)

        #Set output directory + prefix
        sim.use_output_directory(self.logger.data_dir)
        sim.filename_prefix = ['', 'vac'][int(is_vac)]

        #Get fields to output
        if self.save_all:
            fld_names = {'s':['efield_z','hfield_x','hfield_y'], \
                'p':['hfield_z','efield_x','efield_y']}[pol]
        else:
            fld_names = {'s':['efield_z','hfield_y'], \
                'p':['hfield_z','efield_y']}[pol]

        fld_outs = [getattr(mp, f'output_{fn}') for fn in fld_names]

        #Run sim
        sim.run(mp.at_end(mp.synchronized_magnetic(*fld_outs)), \
            until=self.msim.run_time)

        #Output times
        fname = f'{self.logger.data_dir}/{["", "vac-"][int(is_vac)]}times_{pol}'
        sim.output_times(fname)

        #Get metadata
        if get_meta:

            #Get dielectric
            eps = sim.get_array(component=mp.Dielectric)

            #Get coordinates
            x,y,z,w = sim.get_array_metadata()

            #Get run time
            run_time = sim.meep_time()

            #Save metadata
            if semp.zero_rank:

                #Prefix
                pre = f"{self.logger.data_dir}/{['', 'vac-'][int(is_vac)]}"
                pst = f"-{run_time:09.2f}"

                #Save dielectric
                with h5py.File(f'{pre}eps{pst}.h5', 'w') as f:
                    f.create_dataset('dielectric', data=eps)

                #Save coordinates
                with h5py.File(f'{pre}coords{pst}.h5', 'w') as f:
                    f.create_dataset('xx', data=x)
                    f.create_dataset('yy', data=y)
                    f.create_dataset('zz', data=z)

                #Save run time
                np.save(f'{self.logger.data_dir}/time_ext', run_time)

            #Wait
            semp.mpi_barrier()

            #Cleanup
            del eps, x, y, z

        #Reset meep
        sim.reset_meep()

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
