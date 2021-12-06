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
        #Clear taper angle if scallops given
        if len(self.parent.scallop_list) > 0:
            self.parent.taper_angle = 0

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

        #Get initial grid size
        lx0 = 2*(self.padx + self.pmlx) + self.wafer_thick
        if self.is_edge:
            delta_y = self.seam_dark + self.seam_lite
        else:
            delta_y = 2*self.seam_dark + self.gap_width
        ly0 = 2*(self.pady + self.pmly) + delta_y

        #Make sure integer number of pixels
        if (lx0 * self.resolution) % 1 != 0:
            self.padx = np.ceil(lx0*self.resolution)/self.resolution/2 - \
                self.wafer_thick/2 - self.pmlx

        if (ly0 * self.resolution) % 1 != 0:
            self.pady = (np.ceil(ly0*self.resolution)/self.resolution - \
                delta_y - 2*self.pmly) / 2

        #Combine pml and pad
        self.padpmlx = self.padx + self.pmlx
        self.padpmly = self.pady + self.pmly
        self.padpmlz = self.padz + self.pmlz

        ## x ##
        self.lx = 2*self.padpmlx + self.wafer_thick
        self.source_x = -self.lx/2. + self.pmlx
        self.far_pt_x = self.lx/2 - self.pmlx - self.padx/2
        self.non_pml_lx = self.lx - 2*self.pmlx

        ## y ##
        if self.is_edge:
            delta_y = self.seam_dark + self.seam_lite
        else:
            delta_y = 2*self.seam_dark + self.gap_width
        self.ly = 2*self.padpmly + delta_y
        self.edge_y = self.ly/2 - (self.padpmly + self.seam_dark)
        self.non_pml_ly = self.ly - 2*self.pmly

        ## z ##
        self.lz = 0
        self.non_pml_lz = 0

        ## cell ##
        self.cell_size = mp.Vector3(self.lx, self.ly, self.lz)

        #Decay checkpoint
        if self.is_edge:
            self.decay_checkpoint = mp.Vector3(self.wafer_thick/2, self.non_pml_ly/2)
        else:
            self.decay_checkpoint = mp.Vector3(self.wafer_thick/2, -self.gap_width/4)

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
        if zdepth < 10000:
            skcz = zdepth/2
        skn_sze = mp.Vector3(self.skin_thick, sksy, zdepth)
        skn_cen = mp.Vector3(skcx, skcy, skcz)
        skin = mp.Block(material=skn_mat, size=skn_sze, center=skn_cen)
        geometry += [skin]

        #Oxidation layer
        if self.oxide_thick > 0:
            oksy = y0 - y1
            okcx = -(self.wafer_thick/2 + self.skin_thick + self.oxide_thick/2)
            okcy = -y0 + oksy/2.
            okcz = 0
            if not np.isclose(zdepth, mp.inf):
                okcz = zdepth/2
            oxi_sze = mp.Vector3(self.oxide_thick, oksy, zdepth)
            oxi_cen = mp.Vector3(okcx, okcy, okcz)
            oxide = mp.Block(material=oxi_mat, size=oxi_sze, center=oxi_cen)
            geometry += [oxide]

        #Sidewalls
        if self.wall_thick > 0:
            v1 = mp.Vector3(-ex, -y1)
            v0 = mp.Vector3( ex, -y1 - dy)
            dv = mp.Vector3(y=self.wall_thick)
            wal_verts = [v0, v1, v1 - dv, v0 - dv]
            geometry += [mp.Prism(wal_verts, zdepth, material=skn_mat)]

        #Scallops
        if len(self.scallop_list) > 0:
            #Loop through scallop list
            for scallop in self.scallop_list:
                #Get current size and center
                cur_cen, cur_sze = scallop
                scl_sze = mp.Vector3(cur_sze[1], cur_sze[0], zdepth)
                scl_cen = mp.Vector3(-self.wafer_thick/2. + cur_cen[1], -y1+cur_cen[0])
                #Add side walls
                if self.wall_thick > 0:
                    #Side wall ellipsoid
                    geometry += [mp.Ellipsoid(material=skn_mat, size=scl_sze, center=scl_cen)]
                    #Shrink air cylinder's size
                    scl_sze -= mp.Vector3(2*self.wall_thick, 2*self.wall_thick, 0)

                #Add air ellipsoid
                geometry += [mp.Ellipsoid(material=mp.air, size=scl_sze, center=scl_cen)]

            #Clear side walls scallops
            if self.wall_thick > 0:
                #Maximum scallop width + height
                msw = np.max([scl[1][0] for scl in self.scallop_list])
                msh = np.max([scl[1][1] for scl in self.scallop_list])
                #Vertical clear
                p0v = waf_verts[2] - mp.Vector3(self.skin_thick + self.oxide_thick)
                p1v = waf_verts[3] + mp.Vector3(self.wall_thick)
                bsw_vv = [p0v, p0v + mp.Vector3(y=msw), p1v + mp.Vector3(y=msw), p1v]
                geometry += [mp.Prism(bsw_vv, zdepth, material=mp.air)]
                #Horizontal clear
                p0h = waf_verts[3] - mp.Vector3(x=-1/self.parent.resolution,y=msw/2)
                bsw_hv = [p0h, p0h + mp.Vector3(y=msw), p0h + mp.Vector3(msh,msw), \
                    p0h + mp.Vector3(msh)]
                geometry += [mp.Prism(bsw_hv, zdepth, material=mp.air)]

            #Scallop ball
            if self.scallop_ball > 0:
                #Get bottom of first scallop #TODO: add to other scallops
                cur_cen, cur_sze = self.scallop_list[0]
                scl_cen = mp.Vector3(-self.wafer_thick/2. + cur_cen[1] + \
                    cur_sze[1]/2, -y1)
                #Add ball
                geometry += [mp.Cylinder(material=skn_mat, \
                    radius=self.scallop_ball, center=scl_cen, height=zdepth)]

        elif self.scallop_depth > 0:
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

        #Add shaving off of scallops
        if self.shave_angle > 0:
            ex = self.wafer_thick/2.
            dy = self.wafer_thick * (np.tan(np.radians(self.shave_angle)) + \
                np.tan(np.radians(self.taper_angle)))
            #lower (in image) left, upper left, upper right, lower right
            shv_verts = [mp.Vector3( ex, -y1-dy), mp.Vector3(-ex, -y1-1e-3), \
                         mp.Vector3(-ex, -y1),    mp.Vector3( ex, -y1)]
            shave = mp.Prism(shv_verts, float(zdepth), material=mp.air)
            geometry += [shave]

        #Footing
        if self.footing_size is not None:
            fsze = mp.Vector3(x=self.footing_size[0], y=self.footing_size[1], z=zdepth)
            ys0 = -y1 - self.wafer_thick/2*np.tan(np.radians(self.taper_angle))
            fcen = mp.Vector3(x=self.wafer_thick/2, y=ys0)
            geometry += [mp.Ellipsoid(material=mp.air, size=fsze, center=fcen)]

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
