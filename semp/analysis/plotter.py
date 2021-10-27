"""
plotter.py

Author: Anthony Harness
Affiliation: Princeton University
Created on: 08-11-2021
Package: SEMP

Description: Plot results from SEMP simulation
License: Refer to $pkg_home_dir/LICENSE
"""

import semp
import numpy as np
import matplotlib.pyplot as plt;plt.ion()
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import ScalarFormatter

class Plotter(object):

    def __init__(self, alz):
        self.alz = alz

    @property
    def msim(self):
        return self.alz.prop.msim

    @property
    def geo(self):
        return self.alz.prop.msim.geo

############################################
####	Image Plot ####
############################################

    def plot_image(self, data, is_phase=False, title='', vmax=None):

        sx = data.shape[1]/max(data.shape) * 10
        sy = data.shape[0]/max(data.shape) * 10

        fig, axes = plt.subplots(1, figsize=(sx,sy))

        #Get image extent
        extent = [self.alz.yy[0], self.alz.yy[-1], self.alz.xx[-1], self.alz.xx[0]]
        #Adjust for yee lattice
        extent = [x + 0.5/self.msim.resolution for x in extent]

        #Show image
        out = axes.imshow(data, extent=extent, vmax=vmax)

        #Add substrate
        self.draw_substrate(axes)

        #Colorbar
        divider = make_axes_locatable(axes)
        cax = divider.append_axes("right", size="3%", pad=0.05)
        cbar = plt.colorbar(out, cax=cax, format=ScalarFormatter(useMathText=True))
        cbar.ax.yaxis.get_offset_text().set_fontsize(12)
        cbar.set_label(['Amplitude', 'Phase'][int(is_phase)])

        #Labels
        axes.set_xlabel('Y [microns]')
        axes.set_ylabel('X [microns]')
        axes.set_title(title)

        return axes

############################################
############################################

############################################
####	Epsilon Plot ####
############################################

    def plot_epsilon(self, data, with_lines=True):

        fig, axes = plt.subplots(1, figsize=(6,6))

        #Get image extent
        extent = [self.alz.yy[-1], self.alz.yy[0], self.alz.xx[-1], self.alz.xx[0]]
        #Adjust for yee lattice
        extent = [x + 0.5/self.msim.resolution for x in extent]

        #Show image
        out = axes.imshow(data, extent=extent, interpolation='none',
            vmin=0, vmax=data[data < 1e3].max())

        #Add substrate
        if with_lines:
            self.draw_substrate(axes)

        #Labels
        axes.set_xlabel('Y [microns]')
        axes.set_ylabel('X [microns]')

        return axes

############################################
############################################

############################################
####	Slice Plot ####
############################################

    def plot_slice(self, data, is_phase=False):

        fig, axes = plt.subplots(1, figsize=(6,6))

        #Plot
        axes.plot(self.alz.yy, data, '-')
        # axes.plot(self.alz.yy, data[::-1], '--')
        axes.axvline(0, color='k', linestyle=':')

        #Labels
        axes.set_xlabel('Y [microns]')
        axes.set_ylabel(['Amplitude', 'Phase'][int(is_phase)])

        return axes

############################################
############################################

############################################
####	Substrate lines ####
############################################

    def draw_substrate(self, axes):
        #Corners: [top left, top right, bot right, bot left]

        #Widths
        lx = self.geo.lx - 2*self.geo.pmlx
        ly = self.geo.ly - 2*self.geo.pmly
        wx = self.msim.wafer_thick
        wy = -self.geo.edge_y

        #Substrate box
        sub_box = np.array([
            [wy,     -wx/2],
            [-ly/2., -wx/2],
            [-ly/2.,  wx/2],
            [wy,      wx/2]
        ])
        self.overplot_eps_box(sub_box, axes, col='r')

        #Skin box
        if self.msim.skin_thick > 0.:
            skn_box = np.array([
                [wy,     -wx/2 - self.msim.skin_thick],
                [-ly/2,  -wx/2 - self.msim.skin_thick],
                [-ly/2,   wx/2],
                [wy,      wx/2]
            ])
            self.overplot_eps_box(skn_box,axes,col='m')

    def add_box(self, cnrs, axes, col='r'):
        #Loop over corners and plot lines in between
        ncnrs = len(cnrs)
        for i in range(ncnrs):
            #Corner points
            curx = cnrs[:,0][[i,(i+1)%ncnrs]]
            cury = cnrs[:,1][[i,(i+1)%ncnrs]]
            #Shift to edge
            curx += self.geo.edge_y
            #Plot
            axes.plot(curx, cury , color=col)

    def overplot_eps_box(self, cnrs, axes, col='r'):
        #Add box
        self.add_box(cnrs, axes, col=col)

        #Add symmetric box
        if not self.geo.is_edge:
            self.add_box(cnrs*np.array([-1.,1.]), axes, col=col)

############################################
############################################
