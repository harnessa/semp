"""
sommerfeld.py

Author: Anthony Harness
Affiliation: Princeton University
Created on: 12-09-2020
Package: SEMP

Description: Class to calculate Sommerfeld's analytic solution to edge diffraction
License: Refer to $pkg_home_dir/LICENSE
"""

import numpy as np
from scipy.special import fresnel

class Sommerfeld(object):

    def __init__(self, params):
        #Initialize
        self.initialize(params)

    def initialize(self, params, prop=None):
        def_pms = {
            'wave':     0.641,
            'phi0':     np.pi/2,
        }

        #Set default parameters
        for k,v in def_pms.items():
            setattr(self, k, v)

        #Set user parameters
        for k,v in params.items():
            setattr(self, k, v)

        #Derived
        self.kk = 2.*np.pi/self.wave

############################################
####	Analytical Solution ####
############################################

    def G_func(self, s):
        S,C = fresnel(np.sqrt(2./np.pi)*s)
        ans = (1. + 1j)/2 - (C + 1j*S)
        ans *= np.exp(-1j*s**2.)
        del S, C
        return ans

    def get_sommerfeld_solution(self, xx, yy, is_bbek=False):

        #Get coordinates (x=y, y=z)
        rho = np.sqrt(yy**2 + xx**2.)
        if np.isclose(xx, 0):
            phi = (2 - np.heaviside(yy, 1))*np.pi
        else:
            phi = 2.*np.pi + np.arctan2(-xx, -yy)

        #Incident field
        I0 = np.exp(-1j*self.kk*rho*np.cos(phi - self.phi0))

        #Build arguments of calculation
        uu = -np.sqrt(2.*self.kk*rho)*np.cos(0.5*(phi - self.phi0))
        vv = -np.sqrt(2.*self.kk*rho)*np.cos(0.5*(phi + self.phi0))
        pre = (np.exp(-1j*np.pi/4.) / np.sqrt(np.pi) * np.sqrt(np.pi/2)) * \
            np.exp(1j*self.kk*rho)

        #Normalize prefactor by plane wave
        pre /= I0

        #Intermediate calculations
        Umid = pre*self.G_func(uu)
        Gv = pre*self.G_func(vv)
        Dmid = 2j*pre / np.sqrt(np.pi*self.kk*rho)

        #Subtract incident field if Braunbek
        if is_bbek:
            Umid -= np.heaviside(yy, 1)

        #Cleanup
        del uu, vv, pre

        #Get field solution for s,p polarization
        Ez = Umid - Gv                  #s-pol
        Hz = Umid + Gv                  #p-pol

        #Get derivative fields (Meep coords, switched from Born + Wolf; BWx = My, etc.)
        Hy =  np.sin(self.phi0)*Hz + Dmid*np.cos(phi/2)*np.sin(self.phi0/2)     #negative of BW
        Hx = -np.cos(self.phi0)*Ez + Dmid*np.sin(phi/2)*np.sin(self.phi0/2)

        Ey =  np.sin(self.phi0)*Ez + Dmid*np.sin(phi/2)*np.cos(self.phi0/2)
        Ex = -np.cos(self.phi0)*Hz - Dmid*np.cos(phi/2)*np.cos(self.phi0/2)

        #Cleanup
        del Umid, Gv, Dmid, rho, phi, I0

        return Ex, Ey, Ez, Hx, Hy,Hz

############################################
############################################
