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
        self.waves = np.atleast_1d(self.waves)

        #Set pads + pmls
        self.set_pads_pmls()

        #Overwrite Sommerfeld solution
        if self.is_sommerfeld:
            self.apply_sommerfeld()

        #Calculate frequencies and width
        self.freqs = 1./self.waves
        self.fcen0 = self.freqs.mean()
        self.dfreq = max(2.*self.freqs.ptp(), 0.35)

        #Source offset
        self.src_offset = mp.Vector3(y=self.source_offset_y, \
            z=self.source_offset_z) * self.util.m2mu

        #Convert timescales
        if self.prop.save_nt is None:
            self.save_dt = self.n_periods
        else:
            self.save_dt = 1/self.prop.save_nt

        #Force scallop height to be greater than twice depth
        self.scallop_height = max(self.scallop_height, 2*self.scallop_depth)

        #Materials
        self.set_materials()

        #Load Geometry
        if self.sim_geometry == 'corner':
            self.geo = semp.simulation.Geometry_3D(self)
        else:
            self.geo = semp.simulation.Geometry_2D(self)

    def set_pads_pmls(self):
        #Set all components if specified
        for comp in ['x','y','z']:
            if self.pml_all is not None:
                setattr(self, f'pml{comp}', self.pml_all)
            if self.pad_all is not None:
                setattr(self, f'pad{comp}', self.pad_all)

    def apply_sommerfeld(self):
        #Sommerfeld is infinitely thin PEC
        self.wafer_material = 'metal'
        self.skin_material = 'metal'
        self.wafer_epsilon = None
        self.skin_epsilon = None
        for nme in ['skin_thick', 'wall_thick', 'scallop_depth', 'taper_angle']:
            setattr(self, nme, 0)

        #Infinitely thin wafer
        self.wafer_thick = min(0.01, 0.5/self.resolution)

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

            if getattr(self, f'{ob}_epsilon') is not None:
                #Set via epsilon
                eps = getattr(self, f'{ob}_epsilon').real
                #FIXME: choose frequency
                Dcon = 2.*np.pi*self.freqs[0] * getattr(self, f'{ob}_epsilon').imag / \
                    getattr(self, f'{ob}_epsilon').real
                setattr(self, f'{ob}_mat_obj', mp.Medium(epsilon=eps, D_conductivity=Dcon))

            else:
                #Get from library
                setattr(self, f'{ob}_mat_obj', getattr(mat_lib, \
                    getattr(self, f'{ob}_material')))

############################################
############################################

############################################
####	Build Simulation ####
############################################

    def build_sim(self, pol='s', is_vac=False):
        """X is aligned with propagation distance, Y is perpendicular to gap/edge,
           Z is parallel to gap/edge"""

        #K-point
        k_point = mp.Vector3()

        #PML
        pml_layers = self.get_pml_layers(is_vac)

        #Symmetries
        symmetries = self.get_symmetries(is_vac, pol)

        #Computational cell
        cell_size = self.get_cell_size(is_vac)

        #Source
        sources = self.get_source(pol, is_vac)

        #Geometry
        geometry = self.get_geometry(is_vac)

        #Build simulation
        sim = mp.Simulation(split_chunks_evenly=False, force_complex_fields=True,
            ensure_periodicity=False, resolution=self.resolution, Courant=self.courant,
            cell_size=cell_size, boundary_layers=pml_layers, sources=sources,
            geometry=geometry, symmetries=symmetries, k_point=k_point)

        return sim

############################################
############################################

############################################
####	Common Simulation Components ####
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

    def get_symmetries(self, is_vac, pol):
        #Check y symmetry
        if not self.geo.has_y_symm:
            return []

        #No symmetry for vacuum
        if is_vac:
            return []
        else:
            phs = {'s':1, 'p':-1}[pol]
            return [mp.Mirror(mp.Y, phase=phs)]

    def get_source(self, pol, is_vac):
        #Center of source
        src_pt = mp.Vector3(x=self.geo.source_x)

        #Size of source
        src_sze_y = [self.geo.ly, 0.][is_vac]
        src_sze_z = [self.geo.lz, 0.][is_vac]

        #Get source functions
        sim_src, amp_func = self.get_source_function()

        #Get source component
        src_comp = getattr(mp, {'s': 'Ez', 'p': 'Hz'}[pol])

        #Build source
        sources = [mp.Source(sim_src, component=src_comp, center=src_pt, \
            size=mp.Vector3(y=src_sze_y, z=src_sze_z), amp_func=amp_func)]

        return sources

    def get_source_function(self):

        #Get source dependent
        sim_src = mp.GaussianSource(self.fcen0, fwidth=self.dfreq, is_integrated=True)

        #For amp func   #TODO: check amp function with broadband
        kk = 2.*np.pi*self.fcen0
        dist = self.source_distance*self.util.m2mu + self.geo.source_x

        #Amplitude function
        if self.is_diverging:
            def amp_func(pos):
                rr = np.sqrt((pos.y - self.src_offset.y)**2 + \
                    (pos.z - self.src_offset.z)**2 + (pos.x + dist)**2)
                return np.exp(1j*kk*rr)*(dist/rr)

        else:
            if self.src_offset.norm() != 0:
                #Get amplitude function
                amp_func = self.geo.get_src_amp_func(kk)
            else:
                amp_func = None

        return sim_src, amp_func

    def get_geometry(self, is_vac):
        if is_vac:
            return []
        else:
            return self.geo.get_geometry()

############################################
############################################

############################################
####	Misc Functions ####
############################################

    # def adjust_coordinates(self, xx, yy, zz):
    #     #Shift y if edge
    #     if self.sim_geometry == 'edge':
    #         yy += self.gap_width/2
    #
    #     #Shift x to bottom of wafer
    #     # xx -= self.geo.wafer_thick/2
    #
    #     #Shift y by 1 resolution
    #     # yy -= 1/self.resolution
    #
    #     return xx, yy, zz

############################################
############################################
