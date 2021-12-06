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

            #Run sim with wafer
            waf_run_time = self.run_single_to_end(pol, False, None, get_meta)

            #Run sim in vacuum
            vac_run_time = self.run_single_to_end(pol, True, waf_run_time, get_meta)

            #Turn off flag
            get_meta = False

    def run_single_to_end(self, pol, is_vac, run_time, get_meta):

        #Build to simulation
        sim = self.msim.build_sim(pol=pol, is_vac=is_vac)

        #Set output directory + prefix
        sim.use_output_directory(self.logger.data_dir)
        sim.filename_prefix = ['', 'vac'][int(is_vac)]
        dft_name = f'{self.logger.data_dir}/{["", "vac-"][int(is_vac)]}fields_{pol}'

        #Get fields to output
        if self.save_all:
            fld_names = {'s':[mp.Ez, mp.Hx, mp.Hy], 'p':[mp.Hz, mp.Ex, mp.Ey]}[pol]
            far_names = {'s':['Ez', 'Hx', 'Hy'], 'p':['Hz', 'Ex', 'Ey']}[pol]
        else:
            fld_names = {'s':[mp.Ez, mp.Hy], 'p':[mp.Hz, mp.Ey]}[pol]
            far_names = {'s':['Ez', 'Hy'], 'p':['Hz', 'Ey']}[pol]

        #Volume to save (non-PML)
        vlx = self.msim.geo.non_pml_lx
        if is_vac:
            vly, vlz = 0, 0
        else:
            vly = self.msim.geo.non_pml_ly
            vlz = self.msim.geo.non_pml_lz
        vol = mp.Volume(size=mp.Vector3(vlx, vly, vlz))

        #Add DFT fields
        dft_obj = sim.add_dft_fields(fld_names, self.msim.freqs, where=vol)

        #Add near2far
        if self.with_farfield:
            n2f_cen = mp.Vector3(x=self.msim.geo.far_pt_x)
            n2f_sze = mp.Vector3(y=self.msim.geo.non_pml_ly)
            n2f_vol = mp.Near2FarRegion(center=n2f_cen, size=n2f_sze)
            n2f_obj = sim.add_near2far(self.msim.freqs, n2f_vol)

        #Decay check
        dcy_pt = self.msim.geo.decay_checkpoint
        dcy_cn = {'s':mp.Ez, 'p':mp.Hz}[pol]

        #Run sim
        if is_vac:
            #Run vacuum sim to same time as original simulation
            sim.run(until=run_time)

        else:
            #Run wafer sim until decays
            sim.run(until_after_sources=mp.stop_when_fields_decayed( \
                dt=self.msim.decay_dt, pt=dcy_pt, c=dcy_cn, decay_by=self.msim.decay_by))

        #Synchronize magnetic fields
        sim.fields.synchronize_magnetic_fields()

        #Output DFT fields
        sim.output_dft(dft_obj, dft_name)

        #Save simulation time
        sim_time = sim.meep_time()
        sim_timesteps = sim_time / (self.msim.courant / self.msim.resolution)
        if semp.zero_rank:
            tname =  f'{self.logger.data_dir}/{["", "vac-"][int(is_vac)]}simtime_{pol}.h5'
            with h5py.File(tname, 'w') as f:
                f.create_dataset('sim_time', data=sim_time)
                f.create_dataset('sim_timesteps', data=sim_timesteps)
                f.create_dataset('courant', data=self.msim.courant)
                f.create_dataset('resolution', data=self.msim.resolution)

        #Get far fields
        if self.with_farfield:
            far_cen = mp.Vector3(self.msim.farfield_z * self.util.m2mu)
            far_sze = mp.Vector3(y=self.msim.farfield_width * self.util.m2mu)
            ff_res = self.msim.farfield_npts/far_sze.y

            #Get fields
            farfields = sim.get_farfields(n2f_obj, ff_res, center=far_cen, size=far_sze)

            #Save farfields
            if semp.zero_rank:
                n2f_name = f'{self.logger.data_dir}/{["", "vac-"][int(is_vac)]}farfields_{pol}.h5'
                with h5py.File(n2f_name, 'w') as f:
                    f.create_dataset('freqs', data=self.msim.freqs)
                    f.create_dataset('width', data=self.msim.farfield_width)
                    f.create_dataset('distance', data=self.msim.farfield_z)
                    for fn in far_names:
                        f.create_dataset(fn, data=farfields[fn])

            #Cleanup
            del farfields

        #Output times
        tname = f'{self.logger.data_dir}/{["", "vac-"][int(is_vac)]}times_{pol}'
        sim.output_times(tname)

        #Get metadata
        if get_meta:

            #Get dielectric
            eps = sim.get_array(vol=vol, component=mp.Dielectric)
            # eps = np.abs(sim.get_epsilon(self.msim.fcen0, vol=vol))

            #Get coordinates
            x,y,z,w = sim.get_array_metadata(dft_cell=dft_obj)

            #Save metadata
            if semp.zero_rank:

                #Prefix
                pre = f"{self.logger.data_dir}/{['', 'vac-'][int(is_vac)]}"

                #Save dielectric
                with h5py.File(f'{pre}eps.h5', 'w') as f:
                    f.create_dataset('dielectric', data=eps)

                #Save coordinates
                with h5py.File(f'{pre}coords_waves.h5', 'w') as f:
                    f.create_dataset('xx', data=x)
                    f.create_dataset('yy', data=y)
                    f.create_dataset('zz', data=z)
                    f.create_dataset('waves', data=self.msim.waves)

            #Wait
            semp.mpi_barrier()

            #Cleanup
            del eps, x, y, z

        #Reset meep
        sim.reset_meep()

        #Return simulation time
        return sim_time

############################################
############################################

############################################
####	Movie Simulation ####
############################################

    def run_movie(self):
        #Load movie maker
        self.movie_maker = semp.analysis.Movie_Maker(self)

        #Loop over polarizations
        for pol in self.msim.polars:
            #Run sim with wafer
            self.run_single_movie(pol)

    def run_single_movie(self, pol):
        #Is vacuum?
        is_vac = self.msim.sim_geometry == 'vacuum'

        #Build simulation
        sim = self.msim.build_sim(is_vac=is_vac)

        #Set output directory + prefix
        sim.use_output_directory(self.logger.data_dir)
        sim.filename_prefix = ''

        #Get filename
        vac_ext = ['','vac_'][is_vac]
        fld_name = {'s':'ez', 'p':'hz'}[pol]
        src_comp = {'s':mp.Ez, 'p':mp.Hz}[pol]
        filename = f'{vac_ext}movie_{fld_name}'

        #Function to output epsilon
        eps_func = mp.at_beginning(mp.output_epsilon)

        #Function to save field (appended) to h5 file
        if pol == 's':
            fld_func = mp.output_efield_z
        else:
            fld_func = mp.output_hfield_z
        h5_func = mp.to_appended(filename, mp.at_every(self.msim.save_dt, fld_func))

        #Function to save png files
        png_opts = "-Zc dkbluered -a yarg -C $EPS"
        png_func = mp.at_every(self.msim.save_dt, \
            mp.output_png(src_comp, png_opts, rm_h5=True))

        #Output function to save data
        out_funcs = [eps_func, h5_func, png_func]

        #Run sim
        sim.run(*out_funcs, until=self.msim.n_periods)

        #Reset sim
        sim.reset_meep()

        #Create movie
        # self.movie_maker.make_movie(filename)

############################################
############################################
