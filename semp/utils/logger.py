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
        self.log_ext = self.get_log_ext()

    def start_up(self):
        #Create save directory
        if self.writeable:
            self.util.create_directory(self.data_dir)
        #Start
        self.start_time = time.perf_counter()
        self.open_log_file()
        self.save_parameters()
        self.print_starting_message()
        self.log_is_open = True

    def close_up(self):
        #Finish
        if not self.log_is_open:
            return
        self.end_time = time.perf_counter()
        self.print_finishing_message()
        self.close_log_file()
        self.log_is_open = False

############################################
############################################

############################################
#####   Properties #####
############################################

    @property
    def writeable(self):
        return self.prop.do_save and semp.zero_rank and not self.prop.is_analysis

############################################
############################################

############################################
#####  File Functions #####
############################################

    def filename(self,base_name,file_type,ext=None,data_dir=None):
        if data_dir is None:
            data_dir = self.data_dir
        if ext is None:
            ext = self.log_ext
        return f"{data_dir}/{base_name}{ext}.{file_type}"

    def get_log_ext(self, ext=None, pol=None):
        if ext is None:
            ext = self.prop.ext
        if ext != '':
            ext = '__' + ext
        #Add polarization to end of everything
        if pol is None:
            pol = self.prop.msim.polarization
        ext = f'{ext}__{pol}'
        return ext

    def open_log_file(self):
        #Return immediately if not writeable
        if not self.writeable:
            return

        #Open file and register it to close at exit
        self.log_file = open(self.filename('logfile','txt'), 'w')
        atexit.register(self.close_log_file)

    def close_log_file(self):
        #Return immediately if not writeable
        if not self.writeable:
            return

        #Close file if still open
        if not self.log_file.closed:
            self.log_file.close()

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

        #Write to log
        if self.writeable:
            self.log_file.write(new_str)

############################################
############################################

############################################
#####  Messages #####
############################################

    def print_starting_message(self):
        self.write(is_brk=True)
        self.write(txt=f'Starting SEMP run at: {self.prop.session}/{self.prop.ext} ...',is_time=True,n_strs=3)
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
        #Return immediately if not saving or if not zero-rank processor
        if not self.writeable:
            return

        #Save user parameters
        pickle.dump(self.prop.params,      open(self.filename('parameters', 'pck'), 'wb'))
        #Save default parameters
        pickle.dump(semp.utils.def_params, open(self.filename('def_params', 'pck'), 'wb'))

    def save_propagation_data(self, in_data):
        #Let processors catch up
        semp.mpi_barrier()

        #Return immediately if not saving or if not zero-rank processor
        if not self.writeable:
            return

        #Get file name
        fname = self.filename('results', 'h5')

        #Save data
        with h5py.File(fname, 'w') as f:
            #Add metadata
            f.attrs['fld_comp'] = self.prop.msim.src_comp_name
            f.attrs['drv_comp'] = self.prop.msim.drv_comp_name
            #Loop through wafer and vacuum groups
            for grp in ['vac', 'waf']:
                #Create group
                fgrp = f.create_group(grp)
                #Loop through and save data
                for k,v in in_data[grp].items():
                    fgrp.create_dataset(k, data=v)

############################################
############################################

############################################
#####   Loading functions #####
############################################

    def load_parameters(self, alz=None):
        #Load from analyzer
        if alz is not None:
            ext = self.get_log_ext(self, ext=alz.ext, pol=alz.polarization)
            usr_fname = self.filename(self, 'parameters', 'pck', ext=ext, data_dir=alz.data_dir)
            def_fname = self.filename(self, 'def_params', 'pck', ext=ext, data_dir=alz.data_dir)
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

    def load_propagation_data(self):

        #Get file name
        fname = self.filename('results', 'h5')

        #Create data package
        data = {}

        #Save data
        with h5py.File(fname, 'r') as f:
            #Metadata
            data['fld_comp'] = f.attrs['fld_comp']
            data['drv_comp'] = f.attrs['drv_comp']
            #Loop through groups
            for grp in ['vac', 'waf']:
                #Loop through and store data
                for k, v in f[grp].items():
                    data[f'{grp}_{k}'] = v[()]

        return data

    def convert_data_package(self, in_data):
        #Create data package
        data = {}

        #Metadata
        data['fld_comp'] = self.prop.msim.src_comp_name
        data['drv_comp'] = self.prop.msim.drv_comp_name

        #Loop through groups
        for grp in ['vac', 'waf']:
            #Loop through and store data
            for k, v in in_data[grp].items():
                data[f'{grp}_{k}'] = v

        return data

############################################
############################################
