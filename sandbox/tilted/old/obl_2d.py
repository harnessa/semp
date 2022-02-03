import meep as mp
import math
import cmath
import numpy as np
import matplotlib.pyplot as plt;plt.ion()

resolution = 40        # pixels/μm

dpml = 1               # PML thickness
dmat = 10              # length of material

sx = dpml+dmat+dpml
sy = 5

src_comp = [mp.Ez, mp.Hz][0]

cell_size = mp.Vector3(sx,sy,0)
pml_layers = [mp.PML(thickness=dpml,direction=mp.X)]

fcen = 1               # center frequency (wavelength = 1 μm)
df = 0.1               # frequency width

ng = 1.5
glass = mp.Medium(index=ng)

# rotation angle of incident planewave; CCW about Z axis, 0 degrees along +X axis
theta_in = math.radians(22.9)

# k (in source medium) with correct length (plane of incidence: XY)
k = mp.Vector3(math.cos(theta_in),math.sin(theta_in),0).scale(fcen*ng)
if theta_in == 0:
  k = mp.Vector3(0,0,0)

def pw_amp(k,x0):
  def _pw_amp(x):
    return cmath.exp(1j*2*math.pi*k.dot(x+x0))
  return _pw_amp

# src_pt = mp.Vector3(-0.5*sx+dpml+0.3*dmat,0,0)
src_pt = mp.Vector3(-0.5*sx+dpml,0,0)
sources = [mp.Source(mp.GaussianSource(fcen,fwidth=df, is_integrated=True),
                     component=src_comp,
                     center=src_pt,
                     size=mp.Vector3(0,sy,0),
                     amp_func=pw_amp(k,src_pt))]

sim = mp.Simulation(resolution=resolution,
                    cell_size=cell_size,
                    boundary_layers=pml_layers,
                    k_point=k,
                    default_material=glass,
                    sources=sources)

#Add DFT fields
vol = mp.Volume(size=mp.Vector3(sx-2*dpml,sy))
dft_obj = sim.add_dft_fields([src_comp], [fcen], where=vol)

#Run
sim.run(until_after_sources=mp.stop_when_fields_decayed(50, src_comp, \
    mp.Vector3(-0.5*cell_size.x+dpml), 1e-6))

fld = sim.get_dft_array(dft_obj, src_comp, 0)
x,y,z,w = sim.get_array_metadata(dft_cell=dft_obj)


plt.figure()
plt.imshow(np.angle(fld))
plt.figure()
plt.imshow(np.abs(fld))

breakpoint()
