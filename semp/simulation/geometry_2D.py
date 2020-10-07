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

    def __init__(self, parent):
        self.util = semp.utils.Utilities()
        self.parent = parent        # .Meep_Sim
        self.initialize()

############################################
####	Initialization ####
############################################

    def initialize(self):
        #Copy parameters over from parent
        for k in semp.utils.def_params['MEEP_params'].keys():
            setattr(self, k, getattr(self.parent, k))

        #Distance Properties
        self.set_distance_properties()

    def set_distance_properties(self):

        ## cell size ##
        self.lx = self.pmlx + self.padx + self.wafer_thick + self.padx + self.pmlx
        self.ly = self.pmly + self.pady + self.gap_width   + self.pady + self.pmly
        self.lz = 0
        self.cell_size = mp.Vector3(self.lx, self.ly, self.lz)

        ## x ##
        self.source_x = -self.lx/2. + self.pmlx
        self.obs_pt_x = self.wafer_thick/2. + self.obs_distance - 1./self.resolution

        ## y ##
        self.non_pml_sy = self.ly - 2.*self.pmly
        self.obs_pt_y = -self.gap_width/2. + self.obs_distance - 1./self.resolution

############################################
############################################

############################################
####	Dimension Specific ####
############################################

    @property
    def has_y_symm(self):
        #Symmetry is broken by offset source
        src_symm = self.parent.src_offset.norm() == 0
        #Symmetry is broken by edge
        edg_symm = self.sim_geometry != 'edge'
        return src_symm and edg_symm

    def add_pml_layers(self, layers, BLyz):
        layers += [BLyz(thickness=self.pmly, direction=mp.Y)]
        return layers

    def get_src_amp_func(self, kk):
        return lambda pos: np.exp(1j*kk*(pos.y - self.parent.src_offset.y))

############################################
############################################

############################################
####	Geometry ####
############################################

    def get_geometry(self):
        #Build edge on one side
        geometry = self.build_edge()

        #Need to add other side edge if gap and not symmetric
        if self.sim_geometry == 'gap' and not self.has_y_symm:
            #Build edge and flip coordinates
            geometry += self.flip_edge_y(self.build_edge())

        # import pdb;pdb.set_trace()
        return geometry

    def build_edge(self, zdepth=mp.inf):

        #Container
        geometry = []

        #Materials
        waf_mat = self.parent.wafer_mat_obj
        skn_mat = self.parent.skin_mat_obj

        #Wafer
        ex = self.wafer_thick/2.
        ey0 = -self.ly/2
        ey1 = -self.gap_width/2.
        wdy = self.wafer_thick * np.tan(np.radians(self.taper_angle))
        #lower (in image) left, upper left, upper right, lower right
        waf_verts = [mp.Vector3( ex, ey0), mp.Vector3(-ex, ey0), \
                     mp.Vector3(-ex, ey1), mp.Vector3( ex, ey1 - wdy)]
        wafer = mp.Prism(waf_verts, zdepth, material=waf_mat)
        geometry += [wafer]

        #Skin
        if self.skin_thick > 0:
            sksy = (self.ly - self.gap_width)/2.
            skcx = -(self.wafer_thick + self.skin_thick)/2.
            skcy = -self.ly/2. + sksy/2.
            skn_sze = mp.Vector3(self.skin_thick, sksy, zdepth)
            skn_cen = mp.Vector3(skcx, skcy)
            skin = mp.Block(material=skn_mat, size=skn_sze, center=skn_cen)
            geometry += [skin]

        #Sidewalls
        if self.wall_thick > 0:
            v1 = mp.Vector3(-ex, ey1)
            v0 = mp.Vector3( ex, ey1 - ey2)
            dv = mp.Vector3(y=self.wall_thick)
            wal_verts = [v0, v1, v1 - dv, v0 - dv]
            geometry += [mp.Prism(wal_verts, zdepth, material=skn_mat)]

        #Scallops
        if self.scallop_depth > 0:
            #Radius (from height) + Horizontal offset (from depth)
            Rs = self.scallop_height/2.
            ys = Rs - self.scallop_depth
            #Number of scallops
            n_scls = int(np.round(self.wafer_thick / (2*Rs))) + 1
            #Starting vertical center point
            x0 = -self.wafer_thick/2. + Rs
            for i in range(n_scls):
                #Get new horizontal center point (accounting for taper angle)
                y0 = -self.gap_width/2. - (x0 + self.wafer_thick/2) * \
                    np.tan(np.radians(self.taper_angle))
                #Add cylinders
                scl_cen = mp.Vector3(x=x0, y=y0 - ys)
                geometry += [mp.Cylinder(material=mp.air, radius=Rs, \
                    height=zdepth, center=scl_cen)]
                #Step down to next cylinder's position
                x0 += 2*Rs

        return geometry

    def flip_edge_y(self, edge):
        #Loop through shapes
        for cur_shp in edge:

            #Is prism
            if isinstance(cur_shp, mp.geom.Prism):
                cur_shp.center.y *= -1.
                for p in cur_shp.vertices:
                    p.y *= -1.

            #Is Block
            elif isinstance(cur_shp, mp.geom.Block):
                cur_shp.e2 *= -1.
                cur_shp.center.y *= -1.

            #Is Cylinder
            elif isinstance(cur_shp, mp.geom.Cylinder):
                cur_shp.center.y *= -1.

        return edge

############################################
############################################
