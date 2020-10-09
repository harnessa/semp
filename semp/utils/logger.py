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
        self.save_dir = f"{self.prop.save_dir_base}/{self.prop.session_name}"

    def start_up(self):
        #Create save directory
        if self.writeable:
            self.util.create_directory(self.save_dir)
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

    @property
    def ext(self):
        ext = self.prop.save_ext
        if ext != '':
            ext = '__' + ext
        return ext

############################################
############################################

############################################
#####  File Functions #####
############################################

    def filename(self,base_name,file_type,ext=None):
        if ext is None:
            ext = self.ext
        return f"{self.save_dir}/{base_name}{ext}.{file_type}"

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
        self.write(txt=f'Starting SEMP run at: {self.prop.session_name}/{self.prop.save_ext} ...',is_time=True,n_strs=3)
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

        #Save parameters
        pickle.dump(self.prop.params,      open(self.filename('parameters', 'pck'), 'wb'))
        #Save default parameters
        pickle.dump(semp.utils.def_params, open(self.filename('def_params', 'pck'), 'wb'))

############################################
############################################
