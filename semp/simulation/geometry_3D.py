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
        #Dimensions
        self.ndims = 3

        ## cell size ##
        self.lz = self.pmlz + self.padz + self.corner_length + self.pmlz
        self.cell_size = mp.Vector3(self.lx, self.ly, self.lz)

        ## z ##
        self.blk_sze_z = self.pmlz + self.padz
        self.non_pml_sz = self.lz - 2.*self.pmlz

        ## corner ##
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
        eblk = self.build_edge(y0=self.lz/2., y1=self.lz/2. - self.blk_sze_z, \
            zdepth=self.ly)

        #Rotate end block to be aligned with y-axis (rotate moves center to z=0)
        eblk = self.rotate_edge_z2y(eblk)
        #Shift end block to put at end
        eblk = self.shift_edge(eblk, 'z', -(self.lz - self.blk_sze_z)/2.)

        #Build edges
        edge1 = self.build_edge()
        edge2 = self.flip_edge_y(self.build_edge())

        #Shift edges to match end blocks
        edge1 = self.shift_edge(edge1, 'z', -self.lz/2.)
        edge2 = self.shift_edge(edge2, 'z', -self.lz/2.)

        #Add all shapes together (order matters!)
        geometry = edge1 + edge2 + eblk

        #Add broken corner
        if self.with_broken_corner:
            geometry += self.add_broken_corner(eblk, edge1, edge2)

        #Shift by 1 resolution element
        geometry = self.shift_edge(geometry, 'x', -1/self.resolution)
        geometry = self.shift_edge(geometry, 'y', -1/self.resolution)
        geometry = self.shift_edge(geometry, 'z', -1/self.resolution)

        return geometry

    def add_broken_corner(self, eblk, edge1, edge2):
        #Oversize to account for taper angle
        widy = self.ly
        widz = self.lz

        #Build wall with largest edge
        if self.corner_dz > self.corner_dy:

            ### Y Wall ###

            #Get y_wall
            yw_sy = -(edge1[0].vertices[-2].y + self.corner_dy)
            y_wall = self.build_edge(y0=self.ly/2, y1=yw_sy, zdepth=widz)

            #Shift to end block
            yw_dz = -self.lz + self.blk_sze_z + self.corner_dz
            ##Account for taper angle (not if only one wall)
            # yw_dz -= self.wafer_thick * np.arctan(np.radians(self.taper_angle))
            y_wall = self.shift_edge(y_wall, 'z', yw_dz)

            return y_wall

        else:

            ### Z Wall ###

            #Get z_wall and rotate to align along y (resets to 0)
            z_wall = self.rotate_edge_z2y(self.build_edge(y0=self.lz/2, \
                y1=self.lz/2 - widz, zdepth=widy))

            #Shift to edge
            zw_dy = -self.ly/2. - z_wall[0].height/2. + \
                (self.ly/2. - edge2[0].vertices[-2].y) + self.corner_dy
            ##Account for taper angle (not if only one wall)
            # zw_dy -= self.wafer_thick * np.arctan(np.radians(self.taper_angle))
            z_wall = self.shift_edge(z_wall, 'y', zw_dy)

            #Shift to end block
            zw_dz = -self.lz/2 - widz/2 + self.blk_sze_z + self.corner_dz
            z_wall = self.shift_edge(z_wall, 'z', zw_dz)

            return z_wall

############################################
############################################

############################################
####	Shape Manipulation ####
############################################

    def rotate_edge_z2y(self, edge):
        #Get y-centers of each component so that we can realign on wafer block
        ycens = [ee.center.y for ee in edge]

        #Loop through shapes
        for i in range(len(edge)):
            #Current shape
            cur_shp = edge[i]

            #Is Block
            if isinstance(cur_shp, mp.geom.Block):
                cur_shp.e1 = mp.Vector3(1,0,0)
                cur_shp.e2 = mp.Vector3(0,0,1)
                cur_shp.e3 = mp.Vector3(0,1,0)

            #Is Prism
            elif isinstance(cur_shp, mp.geom.Prism):
                for j in range(len(cur_shp.vertices)):
                    cur_shp.vertices[j] = \
                        cur_shp.vertices[j].rotate(mp.Vector3(1,0,0), np.pi/2)
                cur_shp.axis = mp.Vector3(0,1,0)

            #Is Cylinder
            else:
                cur_shp.axis = mp.Vector3(0,1,0)

            #Move center to x=original; y=0; z=around 0, but keeping relative to skin [1]
            cur_shp.center = mp.Vector3(cur_shp.center.x, 0, ycens[i] - ycens[1])

        return edge

############################################
############################################
