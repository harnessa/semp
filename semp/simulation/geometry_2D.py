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

        #Is edge?
        self.is_edge = self.sim_geometry == 'edge'

        #Distance Properties
        self.set_distance_properties()

    def set_distance_properties(self):
        #Dimensions
        self.ndims = 2

        #Combine pml and pad
        self.padpmlx = self.padx + self.pmlx
        self.padpmly = self.pady + self.pmly
        self.padpmlz = self.padz + self.pmlz

        ## x ##
        self.lx = 2*self.padpmlx + self.wafer_thick
        self.source_x = -self.lx/2. + self.pmlx

        ## y ##
        if self.is_edge:
            self.ly = 2*self.padpmly + self.seam_dark + self.seam_lite
        else:
            self.ly = 2*(self.padpmly + self.seam_dark) + self.gap_width
        self.edge_y = self.ly/2 - (self.padpmly + self.seam_dark)

        ## z ##
        self.lz = 0

        ## cell ##
        self.cell_size = mp.Vector3(self.lx, self.ly, self.lz)

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
        edg_symm = not self.is_edge
        return src_symm and edg_symm

    def add_pml_layers(self, layers, BLyz):
        #For Edge, only add to side opposite edge
        if self.is_edge:
            layers += [BLyz(thickness=self.pmly, direction=mp.Y, side=-1)]
        else:
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

        return geometry

    def build_edge(self, y0=None, y1=None, zdepth=mp.inf):
        #Container
        geometry = []

        #Materials
        waf_mat = self.parent.wafer_mat_obj
        skn_mat = self.parent.skin_mat_obj

        #Set start point and width
        if y0 is None:
            y0 = self.ly/2.
        if y1 is None:
            y1 = self.edge_y

        #Wafer
        ex = self.wafer_thick/2.
        dy = self.wafer_thick * np.tan(np.radians(self.taper_angle))
        #lower (in image) left, upper left, upper right, lower right
        waf_verts = [mp.Vector3( ex, -y0), mp.Vector3(-ex, -y0), \
                     mp.Vector3(-ex, -y1), mp.Vector3( ex, -y1 - dy)]
        wafer = mp.Prism(waf_verts, float(zdepth), material=waf_mat)
        geometry += [wafer]

        #Skin
        sksy = y0 - y1
        skcx = -(self.wafer_thick + self.skin_thick)/2.
        skcy = -y0 + sksy/2.
        skcz = 0
        if not np.isclose(zdepth, mp.inf):
            skcz = zdepth/2
        skn_sze = mp.Vector3(self.skin_thick, sksy, zdepth)
        skn_cen = mp.Vector3(skcx, skcy, skcz)
        skin = mp.Block(material=skn_mat, size=skn_sze, center=skn_cen)
        geometry += [skin]

        #Sidewalls
        if self.wall_thick > 0:
            v1 = mp.Vector3(-ex, -y1)
            v0 = mp.Vector3( ex, -y1 - dy)
            dv = mp.Vector3(y=self.wall_thick)
            wal_verts = [v0, v1, v1 - dv, v0 - dv]
            geometry += [mp.Prism(wal_verts, zdepth, material=skn_mat)]

        #Scallops
        if self.scallop_depth > 0:
            #Number of scallops
            n_scls = int(self.wafer_thick / self.scallop_height) + 1
            #Starting vertical center point
            xs0 = -self.wafer_thick/2 + self.scallop_height/2 + self.scallop_start
            #Scallop size
            scl_sze = mp.Vector3(x=self.scallop_height, y=2*self.scallop_depth, z=zdepth)
            #Loop through and add scallops
            for i in range(n_scls):
                #Get new horizontal center point (accounting for taper angle)
                ys0 = -y1 - (xs0 + self.wafer_thick/2) * \
                    np.tan(np.radians(self.taper_angle))
                #Add ellipsoids
                scl_cen = mp.Vector3(x=xs0, y=ys0)
                geometry += [mp.Ellipsoid(material=mp.air, size=scl_sze, center=scl_cen)]
                #Step down to next cylinder's position
                xs0 += self.scallop_height

        return geometry

############################################
############################################

############################################
####	Shape Manipulation ####
############################################

    def shift_edge(self, edge, comp, dshf):
        #Loop through shapes
        for cur_shp in edge:

            #Is prism
            if isinstance(cur_shp, mp.geom.Prism):
                for p in cur_shp.vertices:
                    setattr(p, comp, getattr(p, comp) + dshf)

            #Shift center for all
            setattr(cur_shp.center, comp, getattr(cur_shp.center, comp) + dshf)

        return edge

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
