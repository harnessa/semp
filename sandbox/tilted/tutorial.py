# This example creates an approximate Ez-polarized planewave in vacuum
# propagating at a 45-degree angle, by using a couple of current sources
# with amplitude exp(ikx) corresponding to the desired planewave.
from __future__ import division

import cmath
import math
import meep as mp
import numpy as np
import matplotlib.pyplot as plt;plt.ion()

s = 11  # the size of the computational cell, not including PML
dpml = 1  # thickness of PML layers

sxy = s + 2 * dpml  # cell size, including PML
cell = mp.Vector3(sxy, sxy, 0)

pml_layers = [mp.PML(dpml)]
resolution = 20

rot_angle = np.radians(10)

# pw-amp is a function that returns the amplitude exp(ik(x+x0)) at a
# given point x.  (We need the x0 because current amplitude functions
# in Meep are defined relative to the center of the current source,
# whereas we want a fixed origin.)  Actually, it is a function of k
# and x0 that returns a function of x ...
def pw_amp(k, x0, a):
    def _pw_amp(x):
        return cmath.exp(1j * k.dot(x + x0)) * a
    return _pw_amp

fcen = 1/0.6  # pulse center frequency
df = 0.02  # turn-on bandwidth
kdir = mp.Vector3(np.cos(rot_angle), np.sin(rot_angle))  # direction of k (length is irrelevant)
n = 1 # refractive index of material containing the source
k = kdir.unit().scale(2 * math.pi * fcen * n)  # k with correct length

k_point = mp.Vector3()#fcen*n).rotate(mp.Vector3(z=1), rot_angle)

sources = [
    mp.Source(
        mp.ContinuousSource(fcen, fwidth=df),
        component=mp.Ez,
        center=mp.Vector3(-0.5 * s, 0),
        size=mp.Vector3(0, s),
        amp_func=pw_amp(k, mp.Vector3(x=-0.5 * s), kdir.x)
    ),
    mp.Source(
        mp.ContinuousSource(fcen, fwidth=df),
        component=mp.Ez,
        center=mp.Vector3(0, -0.5 * s),
        size=mp.Vector3(s, 0),
        amp_func=pw_amp(k, mp.Vector3(y=-0.5 * s), kdir.y)
    )
]

sim = mp.Simulation(
    cell_size=cell,
    sources=sources,
    boundary_layers=pml_layers,
    resolution=resolution,
    default_material=mp.Medium(index=n),
    force_complex_fields=True,
    k_point=k_point,
)

t = 400  # run time
# sim.run(mp.at_end(mp.output_efield_z), until=t)
sim.run(until=t)

nonpml_vol = mp.Volume(center=mp.Vector3(), size=mp.Vector3(10,10,0))
sim.plot2D(fields=mp.Ez, output_plane=nonpml_vol)

fld = sim.get_array(mp.Ez)
fig, axes = plt.subplots(1,2)
axes[0].imshow(np.abs(fld).T, origin='lower')
axes[1].imshow(np.angle(fld).T, origin='lower')


breakpoint()
