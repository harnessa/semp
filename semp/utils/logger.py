"""
logger.py

Author:
Affiliation: Princeton University
Created on: 10-06-2020
Package: SEMP

Description: Class to log print outputs during tests or simulations
License: Refer to $pkg_home_dir/LICENSE
"""

import semp
import numpy as np
import atexit
import time
from datetime import datetime
import h5py
import pickle

class Logger(object):

    def __init__(self, prop):
        self.util = semp.utils.util
        self.prop = prop
        #Initialize
        self.initialize()

############################################
#####  Start/Finish #####
############################################

    def initialize(self):
        self.data_dir = f"{self.prop.base_dir}/{self.prop.session}"

    def start_up(self):
        #Create save directory
        self.util.create_directory(self.data_dir)
        #Start
        self.start_time = time.perf_counter()
        self.save_parameters()
        self.print_starting_message()

    def close_up(self):
        #Finish
        self.end_time = time.perf_counter()
        self.print_finishing_message()

############################################
############################################

############################################
#####  File Functions #####
############################################

    def filename(self,base_name,file_type,data_dir=None):
        if data_dir is None:
            data_dir = self.data_dir
        return f"{data_dir}/{base_name}.{file_type}"

############################################
############################################

############################################
#####  Writing Functions #####
############################################

    def write(self, txt='',is_brk=False,n_strs=2,is_time=True,is_err=False,is_skip=True):
        #Return immediately if not zero rank
        if not semp.zero_rank:
            return

        #Build message
        if is_brk:
            new_str = '*'*40 + '\n'
        else:
            if is_time:
                txt += f' ({str(datetime.utcnow())})'
            if is_err:
                txt = self.util.color_string(f'Error! {txt}', self.util.bad_color)
            str_str = '*'*int(n_strs)
            skp = ['','\n'][int(is_skip)]
            new_str = f'{str_str} {txt} {str_str}{skp}'

        #Print
        if self.prop.verbose:
            print(new_str)

############################################
############################################

############################################
#####  Messages #####
############################################

    def print_starting_message(self):
        self.write(is_brk=True)
        self.write(txt=f'Starting SEMP run at: {self.prop.session} ...',is_time=True,n_strs=3)
        self.write(is_brk=True)

    def print_finishing_message(self):
        self.write(is_brk=True)
        self.write(txt='Finished run',n_strs=3)
        self.write(is_brk=True)
        self.write(txt=f'Total Time: {self.end_time - self.start_time:.1f} [s]',\
            is_time=False,n_strs=3)
        self.write(is_brk=True)

############################################
############################################

############################################
#####   Saving functions #####
############################################

    def save_parameters(self):
        #Return immediately if not zero-rank processor
        if not semp.zero_rank:
            return

        #Save user parameters
        pickle.dump(self.prop.params,      open(self.filename('parameters', 'pck'), 'wb'))
        #Save default parameters
        pickle.dump(semp.utils.def_params, open(self.filename('def_params', 'pck'), 'wb'))

############################################
############################################

############################################
#####   Loading functions #####
############################################

    def load_parameters(self, alz=None):
        #Load from analyzer
        if alz is not None:
            usr_fname = self.filename(self, 'parameters', 'pck', data_dir=alz.data_dir)
            def_fname = self.filename(self, 'def_params', 'pck', data_dir=alz.data_dir)
        else:
            usr_fname = self.filename('parameters', 'pck')
            def_fname = self.filename('def_params', 'pck')

        #Load parameters
        usr_pms = pickle.load(open(usr_fname, 'rb'))
        def_pms = pickle.load(open(def_fname, 'rb'))

        #Start with default values
        MEEP_params = def_pms['MEEP_params']
        PROP_params = def_pms['PROP_params']

        #Overwrite with usr params
        for k, v in usr_pms['MEEP_params'].items():
            MEEP_params[k] = v

        for k, v in usr_pms['PROP_params'].items():
            PROP_params[k] = v

        #Join
        params = {'MEEP_params':MEEP_params, 'PROP_params':PROP_params }

        return params

############################################
############################################
