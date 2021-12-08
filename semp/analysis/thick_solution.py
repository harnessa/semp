"""
thick_solution.py

Author: Anthony Harness
Affiliation: Princeton University
Created on: 10-29-2021
Package: SEMP

Description: Class to calculate analytic solution to gap diffraction
License: Refer to $pkg_home_dir/LICENSE
"""

import semp
import numpy as np
from scipy.special import fresnel

class Thick_Solution(object):

    def __init__(self, params):
        self.util = semp.utils.Utilities()
        #Initialize
        self.initialize(params)

    def initialize(self, params, prop=None):
        def_pms = {
            'wave':         0.641,
            'phi0':         np.pi/2.,
            'wafer_thick':  1,
            'gap_width':    5,
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

    def get_thick_solution(self, xx, yy, is_bbek=False):

        import matplotlib.pyplot as plt;plt.ion()   #FIXME: remove
        #Trim to gap
        yy_plt = yy[np.abs(yy) <= self.gap_width/2.]
        yy_gap = yy_plt.copy() + self.gap_width/2

        #Shared calcs
        zz = np.abs(xx) + self.wafer_thick/2.
        kk = 2.*np.pi/self.wave
        H0 = 1
        E0 = 1
        ky = lambda n: np.pi*n/self.gap_width

        beta = lambda n: np.sqrt(kk**2. - ky(n)**2.)
        cos_n = lambda n: np.cos(np.pi*n) - 1.
        exp_n = lambda n: np.exp(1j*beta(n)*zz)
        sqrt_mu_eps = 1.
        om_mu = kk * sqrt_mu_eps
        om_ep = kk / sqrt_mu_eps

        ################

        #TE (s-pol)
        Bn = lambda n: -1j*2./(self.gap_width*beta(n)) * H0 * cos_n(n)
        Ez_func = lambda n:  1j*om_mu / ky(n) * Bn(n) * np.sin(ky(n)*yy_gap) * exp_n(n)
        Hy_func = lambda n: -1j*beta(n) / ky(n) * Bn(n) * np.sin(ky(n)*yy_gap) * exp_n(n)
        Hx_func = lambda n: Bn(n) * np.cos(ky(n)*yy_gap) * exp_n(n)

        # Ez_func2 = lambda n:  4*E0/(n*np.pi) * np.sin(ky(n)*yy_gap) * exp_n(n) * (kk / beta(n))
        # Hy_func2 = lambda n: -4*E0/(n*np.pi) * np.sin(ky(n)*yy_gap) * exp_n(n)
        # Hx_func2 = lambda n: -4j*E0/(self.gap_width*om_mu) * np.cos(ky(n)*yy_gap) * exp_n(n) * (kk / beta(n))

        Ez_func2 = lambda n: 4/(n*np.pi) * (om_mu/beta(n)) * np.sin(ky(n)*yy_gap) * np.exp(1j*beta(n)*zz)
        Hy_func2 = lambda n: 4/(n*np.pi) * np.sin(ky(n)*yy_gap) * np.exp(1j*beta(n)*zz)
        Hx_func2 = lambda n: -4j/(n*np.pi) * (ky(n)**2/beta(n)) * np.cos(ky(n)*yy_gap) * np.exp(1j*beta(n)*zz)

        Nn = int(np.ceil(2.*self.gap_width/self.wave))
        nns = np.arange(1,Nn, 2)            #even n's are zero b/c of cos_n(n)
        nns = nns.reshape(len(nns), 1)

        Ez_tot = Ez_func(nns).sum(0)
        Hx_tot = Hx_func(nns).sum(0)
        Hy_tot = Hy_func(nns).sum(0)


        # Ez_tot2 = Ez_func2(nns).sum(0)
        # Hx_tot2 = Hx_func2(nns).sum(0)
        # Hy_tot2 = Hy_func2(nns).sum(0)
        #
        # fig, axes = plt.subplots(3, 2, figsize=(8,11))
        #
        # axes[0,0].plot(abs(Ez_tot))
        # axes[0,0].plot(abs(Ez_tot2), '--')
        #
        # axes[1,0].plot(abs(Hy_tot))
        # axes[1,0].plot(abs(Hy_tot2), '--')
        #
        # axes[2,0].plot(abs(Hx_tot))
        # axes[2,0].plot(abs(Hx_tot2), '--')
        #
        # axes[0,1].plot(np.angle(Ez_tot))
        # axes[0,1].plot(np.angle(Ez_tot2), '--')
        #
        # axes[1,1].plot(np.angle(Hy_tot))
        # axes[1,1].plot(np.angle(Hy_tot2), '--')
        #
        # axes[2,1].plot(np.angle(Hx_tot))
        # axes[2,1].plot(np.angle(Hx_tot2), '--')
        #
        # breakpoint()

        ################

        #TM (p-pol)
        # arg = lambda n: (n+1.)*np.pi/self.gap_width*yy_plt - n*np.pi/2.
        # An = lambda n: -H0/(om_ep*self.gap_width/2)*(1. + np.cos(np.pi*n))
        # Hz_func = lambda n:  1j * An(n) * om_ep/ky(n) * np.cos(arg(n))
        # Ex_func = lambda n:  1j * An(n) * beta(n)/ky(n) * np.cos(arg(n))
        # Ey_func = lambda n: -1j * An(n) * np.sin(arg(n))
        #
        # Hz_n0 = H0
        # Ex_n0 = Hz_n0 * om_mu / kk
        #
        # Hz_tot = Hz_n0 + Hz_func(nns).sum(0)
        # Ey_tot = Ex_n0 + Ex_func(nns).sum(0)
        # Ex_tot = Ey_func(nns).sum(0)

        An = lambda n: -H0/(om_ep*self.gap_width/2)*(1. + np.cos(np.pi*n))
        Hz_func = lambda n:  1j * An(n) * om_ep/ky(n)   * np.cos(ky(n)*yy_gap) * exp_n(n)
        Ex_func = lambda n:  1j * An(n) * beta(n)/ky(n) * np.cos(ky(n)*yy_gap) * exp_n(n)
        Ey_func = lambda n: An(n) * np.sin(ky(n)*yy_gap) * exp_n(n)

        Hz_n0 = H0
        Ex_n0 = Hz_n0 * om_mu / kk

        Hz_tot = Hz_n0 + Hz_func(nns).sum(0)
        Ey_tot = Ex_n0 + Ex_func(nns).sum(0)
        Ex_tot = Ey_func(nns).sum(0)

        ################

        #Pad with zeros
        def pad_func(inarr):
            new = np.zeros_like(yy) + 0j
            new[np.abs(yy) <= self.gap_width/2.] = inarr.copy()
            return new

        Ex = pad_func(Ex_tot)
        Ey = pad_func(Ey_tot)
        Ez = pad_func(Ez_tot)
        Hx = pad_func(Hx_tot)
        Hy = pad_func(Hy_tot)
        Hz = pad_func(Hz_tot)

        #Normalize by plane wave at that point
        pw = -np.exp(1j*kk*zz)      #FIXME why negative?
        Ex /= pw
        Ey /= pw
        Ez /= pw
        Hx /= pw
        Hy /= pw
        Hz /= pw

        #Hy needs -1 multiply
        Hy *= np.exp(1j*np.pi)

        return Ex, Ey, Ez, Hx, Hy, Hz

############################################
############################################
