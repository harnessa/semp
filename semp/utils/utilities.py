"""
utilities.py

Author: Anthony Harness
Affiliation: Princeton University
Created on: 10-06-2020
Package: SEMP

Description: Utility functions to be used by SEMP
License: Refer to $pkg_home_dir/LICENSE
"""

import numpy as np
import copy
import os

class Utilities(object):

    def __init__(self):
        ### Text Colors ###
        self.colors = ['r','g','b','y','m','c','k']
        self.bad_color = "\x1b[0;30;41m"
        self.good_color = "\x1b[0;30;42m"
        self.off_color = '\x1b[0m'
        self.neutral_color = "\x1b[0;33;40m"

    def deepcopy(self,in_obj):
        return copy.deepcopy(in_obj)

    def create_directory(self,directory):
        if not os.path.exists(directory):
            os.makedirs(directory)
        return directory

    def color_string(self,string,color):
        return color + string + self.off_color

############################################
#####   Parameter Functions    #####
############################################

    def mix_usr_def_params(self, usr_pms, def_pms):
        usr_keys = usr_pms.keys()
        for k,v in def_pms.items():
            if k not in usr_keys:
                usr_pms[k] = v
        return usr_pms

    def set_default_params(self, parent, params, def_pms):
        bad_str = self.color_string('!*!', self.bad_color)

        #Set default parameters
        for k,v in def_pms.items():
            setattr(parent, k, v)

        #Overwrite with user-specified parameters
        def_keys = def_pms.keys()
        for k,v in params.items():
            #Error message if unknown parameter supplied
            if k not in def_keys:
                print('\n' + bad_str + ' New Parameter not in Defaults: %s ' \
                    %self.color_string(k, self.neutral_color) + bad_str + '\n')
                import sys; sys.exit(0)

            setattr(parent, k, v)

############################################
############################################

util = Utilities()
