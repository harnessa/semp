"""
geometry_3D.py

Author: Anthony Harness
Affiliation: Princeton University
Created on: 10-07-2020
Package: SEMP

Description: Class to build 3D geometry
License: Refer to $pkg_home_dir/LICENSE
"""

import semp
import numpy as np
import meep as mp
from semp.simulation import Geometry_2D

class Geometry_3D(Geometry_2D):

    def __init__(self, parent, params):
        Geometry_2D.__init__(self, parent, params)
        self.initialize_3D()

############################################
####	Initialization ####
############################################

    def initialize_3D(self):
        #Distance Properties
        self.set_distance_properties_3D()

    def set_distance_properties_3D(self):

        ## cell size ##
        self.lz = self.pmlz + self.padz + self.corner_length + self.pmlz
        self.cell_size = mp.Vector3(self.lx, self.ly, self.lz)

        ## z ##
        self.non_pml_sz = self.lz - 2.*self.pmlz
        self.hole_sze_z = self.corner_length + self.pmlz - self.corner_dz - self.wall_thick
        self.hole_cen_z = self.lz/2. - self.hole_sze_z/2
        self.with_broken_corner = self.corner_dy > 0 and self.corner_dz > 0

############################################
############################################

############################################
####	Geometry Specific ####
############################################

    def add_pml_layers(self, layers, BLyz):
        layers.append(BLyz(thickness=self.pmly, direction=mp.Y))
        layers.append(BLyz(thickness=self.pmlz, direction=mp.Z))
        return layers

    def get_symmetries(self):
        #Symmetry is broken by broken corner
        if self.with_broken_corner:
            return []
        else:
            return Geometry_2D.get_symmetries(self)

############################################
############################################
