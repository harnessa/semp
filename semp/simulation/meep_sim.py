"""
meep_sim.py

Author: Anthony Harness
Affiliation: Princeton University
Created on: 10-06-2020
Package: SEMP

Description: Class to hold and run Meep simulation
License: Refer to $pkg_home_dir/LICENSE
"""

import semp
import numpy as np
import meep as mp
import meep.materials as mat_lib

class Meep_Sim(object):

    def __init__(self, prop, params):
        self.util = semp.utils.Utilities()
        self.prop = prop
        #Initialize
        self.initialize(params)

############################################
####	Initialization ####
############################################

    def initialize(self, params):
        #Set meep parameters
        self.util.set_default_params(self, params, semp.utils.def_params['MEEP_params'])

        #Force datatypes
        self.resolution = int(self.resolution)

        #Calculate central frequency and width
        self.fcen = 1./self.wave
        self.dfreq = self.fwid*self.fcen

        #Polarization
        self.src_comp_name = {'s': 'Ez', 'p': 'Hz'}[self.polarization]
        self.src_comp = getattr(mp, self.src_comp_name)

        # #Convert timescales
        # self.run_time = self.n_periods*self.wave
        # if self.save_nt is None:
        #     self.save_dt = self.run_time
        # else:
        #     self.save_dt = self.wave/self.save_nt

        #Materials
        self.set_materials()

        #Load Geometry
        if self.sim_geometry == 'corner':
            self.geo = semp.simulation.Geometry_3D(self, params)
        else:
            self.geo = semp.simulation.Geometry_2D(self, params)
        print(self.geo.get_symmetries())
        import pdb;pdb.set_trace()

############################################
############################################

############################################
####	Materials ####
############################################

    def set_materials(self):

        #Load materials library
        mat_lib.metal = mp.metal

        #Set material objects
        for ob in ['wafer', 'skin']:

            if getattr(self, '%s_epsilon'%ob) is not None:
                #Set via epsilon
                eps = getattr(self, '%s_epsilon'%ob).real
                Dcon = 2.*np.pi*self.fcen * getattr(self, '%s_epsilon'%ob).imag / \
                    self.wafer_epsilon.real
                setattr(self, '%s_mat_obj'%ob, mp.Medium(epsilon=eps, D_conductivity=Dcon))

            else:
                #Get from library
                setattr(self, '%s_mat_obj'%ob, getattr(mat_lib, \
                    getattr(self,'%s_material'%ob)))

############################################
############################################

############################################
####	Build Simulation ####
############################################

    def build_sim(self, is_vac=False):
        """X is aligned with propagation distance, Y is perpendicular to gap/edge,
           Z is parallel to gap/edge"""

        #K-point
        k_point = mp.Vector3()

        #PML
        pml_layers = self.get_pml_layers(is_vac)

        #Symmetries
        symmetries = self.get_symmetries(is_vac)

        #Computational cell
        cell_size = self.get_cell_size(is_vac)

        #Source
        sources = self.get_source(is_vac)

        #Geometry
        geometry = self.get_geometry(is_vac)

        #Build simulation
        sim = mp.Simulation(split_chunks_evenly=True, force_complex_fields=True,
            ensure_periodicity=False, resolution=self.resolution, Courant=self.courant,
            cell_size=cell_size, boundary_layers=pml_layers, sources=sources,
            geometry=geometry, symmetries=symmetries, k_point=k_point)

        return sim

############################################
############################################

############################################
####	Common Wrappers ####
############################################

    def get_cell_size(self, is_vac):
        if is_vac:
            return mp.Vector3(self.geo.lx)
        else:
            return self.geo.cell_size

    def get_pml_layers(self, is_vac):
        #Always use PML for X direction (don't have to worry about entering metal)
        layers = [mp.PML(thickness=self.pmlx, direction=mp.X)]
        #Add other BL if not vacuum
        if not is_vac:
            #Use absorber or PML as boundary layer
            if self.use_absorber:
                BLyz = mp.Absorber
            else:
                BLyz = mp.PML
            layers = self.geo.add_pml_layers(layers, BLyz)
        return layers

    def get_symmetries(self, is_vac):
        #Symmetry is broken by offset source
        if self.source_offset != 0:
            return []

        #No symmetry for vacuum
        if is_vac:
            return []
        else:
            return self.geo.get_symmetries()

    def get_source(self, is_vac):
        #Center of source
        src_pt = mp.Vector3(x=self.source_x)

        #Size of source
        src_sze_y = [self.geo.ly, 0.][is_vac]
        src_sze_z = [self.geo.lz, 0.][is_vac]

        #Get source functions
        sim_src, amp_func = self.get_source_function()

        #Build source   #TODO: add gaussian beam source option
        sources = [mp.Source(sim_src, component=self.src_comp, center=src_pt, \
            size=mp.Vector3(y=src_sze_y, z=src_sze_z), amp_func=amp_func)]

        return sources

    def get_source_function(self):

        #Get source dependent
        sim_src = mp.ContinuousSource(self.fcen, fwidth=self.df, is_integrated=True)

        #For amp func
        kk = 2.*np.pi*self.fcen
        y0 = self.source_offset*self.util.m2mu
        zz = self.source_distance*self.util.m2mu

        #Amplitude function     #TODO: fix for corner
        if self.is_diverging:
            def amp_func(pos):
                rr = np.sqrt((pos.y - y0)**2 + (pos.x + zz)**2)
                return np.exp(1j*kk*rr)*(zz/rr)

        else:
            if self.source_y0 != 0:
                amp_func = lambda pos: np.exp(1j*kk*(pos.y - y0))
            else:
                amp_func = None

        return sim_src, amp_func

############################################
############################################
