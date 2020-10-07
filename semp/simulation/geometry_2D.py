"""
geometry_2D.py

Author: Anthony Harness
Affiliation: Princeton University
Created on: 10-07-2020
Package: SEMP

Description: Class to build 2D geometry
License: Refer to $pkg_home_dir/LICENSE
"""

import semp
import numpy as np
import meep as mp

class Geometry_2D(object):

    def __init__(self, parent, params):
        self.util = semp.utils.Utilities()
        self.parent = parent        # .Meep_Sim
        self.initialize(params)
        # import pdb;pdb.set_trace()

############################################
####	Initialization ####
############################################

    def initialize(self, params):
        self.util.set_default_params(self, params, semp.utils.def_params['MEEP_params'])

        #Distance Properties
        self.set_distance_properties()

    def set_distance_properties(self):

        ## cell size ##
        self.lx = self.pmlx + self.padx + self.wafer_thick + self.padx + self.pmlx
        self.ly = self.pmly + self.pady + self.gap_width   + self.pady + self.pmly
        self.lz = 0
        self.cell_size = mp.Vector3(self.lx, self.ly)

        ## x ##
        self.source_x = -self.lx/2. + self.pmlx
        self.obs_pt_x = self.wafer_thick/2. + self.obs_distance - 1./self.resolution

        ## y ##
        self.non_pml_sy = self.ly - 2.*self.pmly
        self.shadow_y = self.gap_width/2. - self.wall_thick
        self.obs_pt_y = -self.shadow_y + self.obs_distance - 1./self.resolution

############################################
############################################

############################################
####	Geometry Specific ####
############################################

    def add_pml_layers(self, layers, BLyz):
        layers.append(BLyz(thickness=self.pmly, direction=mp.Y))
        return layers

    def get_symmetries(self):
        phs = {'s':1, 'p':-1}[self.polarization]
        # return [mp.Mirror(mp.Y, phase=phs)]
        return phs

############################################
############################################
