"""
movie_maker.py

Author: Anthony Harness
Affiliation: Princeton University
Created on: 10-09-2020
Package: SEMP

Description: Class to make movie out of simulation data
License: Refer to $pkg_home_dir/LICENSE
"""

import semp
import numpy as np
import subprocess
import h5py

class Movie_Maker(object):

    def __init__(self, parent):
        self.parent = parent        #Propagator
        self.tmp_dir = f"{semp.int_data_dir}/tmp"

############################################
####	Main Script ####
############################################

    def make_movie(self, file_end):
        #Wait for processors to catch up
        semp.MPI.COMM_WORLD.Barrier()

        #Return if not zero rank
        if not (self.parent.do_save and semp.zero_rank):
            return

        #Create movie
        self.run_make_movie(file_end)

        #Delete images
        self.run_clean_images()

############################################
############################################

############################################
####	Creation Functions ####
############################################

    def make_images_from_h5(self, filename, w_eps=False):
        #Get size of time slice
        with h5py.File(f'{filename}.h5', 'r') as f:
            n_imgs = f[next(iter(f.keys()))].shape[-1] - 1

        #With epsilon
        eps = ['', '-C $EPS'][w_eps]

        #Build command
        cmd = f'h5topng -t 0:{n_imgs} -R -Zc dkbluered -a yarg {eps} {filename}.h5'

        #Run command and return status
        return self.run_shell_command(cmd)

    def run_make_movie(self, file_end):
        #File names
        png_name = f'{self.parent.logger.save_dir}/{self.parent.meep_sim.src_comp_name.lower()}'
        gif_name = f'{self.parent.logger.save_dir}/{file_end}'

        #Build command
        cmd = f'convert {png_name}*.png {gif_name}.gif'

        #Run command and return status
        return self.run_shell_command(cmd)

    def run_clean_images(self):
        #File names
        png_name = f'{self.parent.logger.save_dir}/{self.parent.meep_sim.src_comp_name.lower()}'

        #Build command
        cmd = f'rm {png_name}*.png'

        #Run command and return status
        return self.run_shell_command(cmd)

############################################
############################################

############################################
####	Common Functions ####
############################################

    def run_shell_command(self, cmd):
        #Run command
        out = subprocess.run(cmd, shell=True, capture_output=True)

        #Check returncode
        if out.returncode != 0:
            cmd_func = cmd.split(' ')[0]
            print(f'\nWARNING! Subprocess: {cmd_func} FAILED with error: {out}!\n')
            return False

        return True

############################################
############################################
