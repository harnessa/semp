import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

import meep as mp
from meep.materials import SiO2, cSi

resolution = 50  # pixels/um

dpml = 1.0
pml_layers = [mp.PML(thickness=dpml)]

r = 1.0     # radius of cylinder
dair = 2.0  # air padding thickness

s = 2*(dpml+dair+r)
cell_size = mp.Vector3(s,s)

wvl = 0.6
fcen = 1/wvl

mat_name = ['SiO2','cSi'][1]

mat = {'SiO2':SiO2,'cSi':cSi}[mat_name]

# is_integrated=True necessary for any planewave source extending into PML
sources = [mp.Source(mp.GaussianSource(fcen,fwidth=0.1*fcen,is_integrated=True),
                     center=mp.Vector3(-0.5*s+dpml),
                     size=mp.Vector3(0,s),
                     component=mp.Ez)]

symmetries = [mp.Mirror(mp.Y)]

# geometry = [mp.Cylinder(material=mat,
#                         center=mp.Vector3(),
#                         radius=r,
#                         height=mp.inf)]

geometry = []#mp.Block(material=mat,
                        # center=mp.Vector3(),
                        # size=mp.Vector3(2*r,mp.inf,mp.inf))]

sim = mp.Simulation(resolution=resolution,
                    cell_size=cell_size,
                    boundary_layers=pml_layers,
                    sources=sources,
                    k_point=mp.Vector3(),
                    symmetries=symmetries,
                    geometry=geometry)

dft_fields = sim.add_dft_fields([mp.Dz,mp.Ez],
                                fcen,0,1,
                                center=mp.Vector3(),
                                size=mp.Vector3(2*r,2*r),
                                yee_grid=True)

# closed box surrounding cylinder for computing total incoming flux
# flux_box = sim.add_flux(fcen, 0, 1,
#                         mp.FluxRegion(center=mp.Vector3(x=-r),size=mp.Vector3(0,2*r),weight=+1),
#                         mp.FluxRegion(center=mp.Vector3(x=+r),size=mp.Vector3(0,2*r),weight=-1),
#                         mp.FluxRegion(center=mp.Vector3(y=+r),size=mp.Vector3(2*r,0),weight=-1),
#                         mp.FluxRegion(center=mp.Vector3(y=-r),size=mp.Vector3(2*r,0),weight=+1))

box1 = sim.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(x=-r),size=mp.Vector3(0,2*r),weight=+1))
box2 = sim.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(x=+r),size=mp.Vector3(0,2*r),weight=-1))

sim.run(until_after_sources=100)

box1_data = sim.get_flux_data(box1)
box2_data = sim.get_flux_data(box2)
box1_flux0 = mp.get_fluxes(box1)
box2_flux0 = mp.get_fluxes(box2)

sim.reset_meep()

geometry = [mp.Block(material=mat,
                        center=mp.Vector3(),
                        size=mp.Vector3(2*r,mp.inf,mp.inf))]

sim = mp.Simulation(resolution=resolution,
                    cell_size=cell_size,
                    boundary_layers=pml_layers,
                    sources=sources,
                    k_point=mp.Vector3(),
                    symmetries=symmetries,
                    geometry=geometry)

box1 = sim.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(x=-r),size=mp.Vector3(0,2*r),weight=+1))
box2 = sim.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(x=+r),size=mp.Vector3(0,2*r),weight=-1))

sim.load_minus_flux_data(box1, box1_data)
sim.load_minus_flux_data(box2, box2_data)

sim.run(until_after_sources=100)

box1_flux = mp.get_fluxes(box1)
box2_flux = mp.get_fluxes(box2)

scatt_flux = np.asarray(box1_flux)-np.asarray(box2_flux)
intensity = np.asarray(box1_flux0)/(2*r)**2
scatt_cross_section = np.divide(scatt_flux,intensity)
scatt_eff_meep = scatt_cross_section*-1/(np.pi*r**2)

alpha = np.log(scatt_cross_section)/(2*r)   #aspnes and studna: 0.5738 /micron at 0.6 um
breakpoint()
Dz = sim.get_dft_array(dft_fields,mp.Dz,0)
Ez = sim.get_dft_array(dft_fields,mp.Ez,0)
absorbed_power_density = 2*np.pi*fcen * np.imag(np.conj(Ez)*Dz)

dxy = 1/resolution**2
absorbed_power = np.sum(absorbed_power_density)*dxy
absorbed_flux = mp.get_fluxes(flux_box)[0]
err = abs(absorbed_power-absorbed_flux)/absorbed_flux
print("flux:, {} (dft_fields), {} (dft_flux), {} (error)".format(absorbed_power,absorbed_flux,err))

plt.figure()
sim.plot2D()
plt.savefig(f'power_density_cell_{mat_name}.png',dpi=150,bbox_inches='tight')

plt.figure()
x = np.linspace(-r,r,Dz.shape[0])
y = np.linspace(-r,r,Dz.shape[1])
plt.pcolormesh(x,
               y,
               np.transpose(absorbed_power_density),
               cmap='inferno_r',
               shading='gouraud',
               vmin=0,
               vmax=np.amax(absorbed_power_density))
plt.xlabel("x (μm)")
plt.xticks(np.linspace(-r,r,5))
plt.ylabel("y (μm)")
plt.yticks(np.linspace(-r,r,5))
plt.gca().set_aspect('equal')
plt.title("absorbed power density" + "\n" +"SiO2 Labs(λ={} μm) = {:.2f} μm".format(wvl,wvl/np.imag(np.sqrt(mat.epsilon(fcen)[0][0]))))
plt.colorbar()
plt.savefig(f'power_density_map_{mat_name}.png',dpi=150,bbox_inches='tight')
breakpoint()
