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

    def __init__(self, parent):
        Geometry_2D.__init__(self, parent)
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
        self.blk_sze_z = self.pmlz + self.padz
        self.with_broken_corner = self.corner_dy > 0 and self.corner_dz > 0

############################################
############################################

############################################
####	Dimension Specific ####
############################################

    @property
    def has_y_symm(self):
        #Symmetry is broken by offset source
        src_symm = self.parent.src_offset.norm() == 0
        #Symmetry is broken by broken corner
        cnr_symm = not self.with_broken_corner
        return src_symm and cnr_symm

    def add_pml_layers(self, layers, BLyz):
        layers = Geometry_2D.add_pml_layers(self, layers, BLyz)
        layers.append(BLyz(thickness=self.pmlz, direction=mp.Z))
        return layers

    def get_src_amp_func(self, kk):
        #FIXME: not sure if this works!
        #Get angle of incidence
        inc_ang = np.arctan2(self.parent.source_offset_y, self.parent.source_offset_z)
        #Get propagation vector (in YZ plane)
        kr = mp.Vector3(kk).rotate(mp.Vector3(x=1), inc_ang)
        #Get offset
        off = mp.Vector3(self.source_x, self.parent.src_offset.y, self.parent.src_offset.z)
        #Build amplitude function
        amp_func = lambda pos: np.exp(1j*kr.dot(pos - off))
        return amp_func

############################################
############################################

############################################
####	Geometry ####
############################################

    def get_geometry(self):

        #Build end block
        end_blk = self.build_end_block()

        #Build edges
        edge1 = self.build_edge()
        edge2 = self.flip_edge_y(self.build_edge())

        #Shift edges to match end block
        dz = self.lz/2.
        edge1 = self.shift_edge_z(edge1, -dz)
        edge2 = self.shift_edge_z(edge2, -dz)

        #Add all shapes together (order matters!)
        geometry = edge1 + edge2 + end_blk

        # import pdb;pdb.set_trace()
        return geometry

    def build_end_block(self):

        #Container
        geometry = []

        #Materials
        waf_mat = self.parent.wafer_mat_obj
        skn_mat = self.parent.skin_mat_obj

        #Wafer
        ex = self.wafer_thick/2.
        ey = -self.ly/2.
        ez0 = -self.lz/2
        ez1 = ez0 + self.blk_sze_z
        wdz = self.wafer_thick * np.tan(np.radians(self.taper_angle))
        #lower (in image) left, upper left, upper right, lower right
        waf_verts = [mp.Vector3(x= ex, y=ey, z=ez0), mp.Vector3(x=-ex, y=ey, z=ez0), \
                     mp.Vector3(x=-ex, y=ey, z=ez1), mp.Vector3(x= ex, y=ey, z=ez1 - wdz)]
        wafer = mp.Prism(waf_verts, self.ly, axis=mp.Vector3(0,1,0), material=waf_mat)
        geometry += [wafer]

        #Skin
        if self.skin_thick > 0:
            skcx = -(self.wafer_thick + self.skin_thick)/2.
            skcz = (-self.lz + self.blk_sze_z)/2.
            skn_sze = mp.Vector3(self.skin_thick, self.ly, self.blk_sze_z)
            skn_cen = mp.Vector3(x=skcx, z=skcz)
            skin = mp.Block(material=skn_mat, size=skn_sze, center=skn_cen)
            geometry += [skin]

        return geometry

    def shift_edge_z(self, edge, dz):
        #Loop through shapes
        for cur_shp in edge:

            #Is prism
            if isinstance(cur_shp, mp.geom.Prism):
                cur_shp.center.z += dz
                for p in cur_shp.vertices:
                    p.z += dz

            #Is Block
            elif isinstance(cur_shp, mp.geom.Block):
                cur_shp.center.z += dz

            #Is Cylinder
            elif isinstance(cur_shp, mp.geom.Cylinder):
                cur_shp.center.z += dz

        return edge

############################################
############################################
