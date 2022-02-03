import numpy as np
import meep as mp
import meep.materials as mat_lib
import h5py
import os

############################################
####	User Input ####
############################################

#Numerics
resolution = 30
dpml = 4
dpad = 4
seam_dark = 6

#Wafer
gap_width = 11
skin_thick = 0.3
device_thick = 2.25
support_thick = 50
lip_width = 10
wafer_material = 'cSi'
skin_material = 'Al'

session = 'test'

#Directory
save_dir = f'/home/aharness/Research/Optics_Modeling/Semp_Results/full_wafer/{session}'
# save_dir = f'/scratch/network/aharness/Semp_Results/full_wafer/{session}'

if mp.am_master():
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

############################################
####	Other params ####
############################################

#Decay
decay_dt = 30         # Time after decay to
decay_by = 1e-6       # Decay amount
waves = np.array([0.641, 0.660, 0.699, 0.725])

#Derived
freqs = 1./waves
fcen0 = freqs.mean()
dfreq = max(2.*freqs.ptp(), 0.35)

mat_lib.metal = mp.metal
waf_mat = getattr(mat_lib, wafer_material)
skn_mat = getattr(mat_lib, skin_material)

lx_np = 2*dpad + support_thick + device_thick
ly_np = 2*dpad + 2*seam_dark + 2*lip_width + gap_width
cell_size = mp.Vector3(dpml*2 + lx_np, dpml*2 + ly_np)

############################################
####	Simulation Function ####
############################################

def run_sim(pol, is_vac, get_meta):

    #Print
    if mp.am_master():
        print(f'\nRunning: {pol}, {is_vac}\n')

    #Run specific
    src_comp = getattr(mp, {'s': 'Ez', 'p': 'Ey'}[pol])

    #Source
    sim_src = mp.GaussianSource(fcen0, fwidth=dfreq, is_integrated=True)
    sources = [mp.Source(sim_src, component=src_comp,
                center=mp.Vector3(x=-0.5*lx_np), \
                size=mp.Vector3(y=cell_size.y))]

    #PML
    pml_layers = [mp.PML(thickness=dpml, direction=mp.X), \
        mp.Absorber(thickness=dpml, direction=mp.Y)]

    #Symmetry
    symmetries = [mp.Mirror(mp.Y, phase={'s':1, 'p':-1}[pol])]

    #Geometry
    geometry = []
    if not is_vac:

        #Total height
        xTot = support_thick + device_thick

        #Support wafer
        sup_sz = mp.Vector3(support_thick, dpml + dpad + seam_dark, mp.inf)
        sup_cc = mp.Vector3(xTot/2 - support_thick/2, -0.5*cell_size.y + 0.5*sup_sz.y)
        support_wafer_L = mp.Block(material=waf_mat, size=sup_sz, center=sup_cc)
        support_wafer_R = mp.Block(material=waf_mat, size=sup_sz, center=sup_cc)
        support_wafer_R.center.y *= -1
        geometry += [support_wafer_L, support_wafer_R]

        #Device wafer
        dev_sz = mp.Vector3(device_thick, sup_sz.y + lip_width, mp.inf)
        dev_cc = mp.Vector3(sup_cc.x - support_thick/2 - device_thick/2, \
            -0.5*cell_size.y + 0.5*dev_sz.y)
        device_wafer_L = mp.Block(material=waf_mat, size=dev_sz, center=dev_cc)
        device_wafer_R = mp.Block(material=waf_mat, size=dev_sz, center=dev_cc)
        device_wafer_R.center.y *= -1
        geometry += [device_wafer_L, device_wafer_R]

        #Bottom of wafers
        device_bottom = device_wafer_L.center.x + device_wafer_L.size.x/2
        support_bottom = support_wafer_L.center.x + support_wafer_L.size.x/2

        #Skin
        skn_sz = mp.Vector3(skin_thick, dev_sz.y, mp.inf)
        skn_cc = mp.Vector3(dev_cc.x - device_thick/2 - skin_thick/2, dev_cc.y)
        skin_L = mp.Block(material=skn_mat ,size=skn_sz, center=skn_cc)
        skin_R = mp.Block(material=skn_mat ,size=skn_sz, center=skn_cc)
        skin_R.center.y *= -1
        geometry += [skin_L, skin_R]

    breakpoint()
    #Build simulation
    sim = mp.Simulation(force_complex_fields=True,
        resolution=resolution, ensure_periodicity=False, eps_averaging=False,
        cell_size=cell_size, boundary_layers=pml_layers, sources=sources,
        geometry=geometry, k_point=mp.Vector3(), symmetries=symmetries)

    ############################################
    ############################################

    #Set output directory + prefix
    sim.use_output_directory(save_dir)
    sim.filename_prefix = ['', 'vac'][int(is_vac)]
    dft_name = f'{save_dir}/{["", "vac-"][int(is_vac)]}fields_{pol}'

    #Get fields to output
    fld_names = {'s':[mp.Ez], 'p':[mp.Ey]}[pol]

    #Add DFT fields
    vol = mp.Volume(size=mp.Vector3(lx_np, ly_np))
    dft_obj = sim.add_dft_fields(fld_names, freqs, where=vol)

    #Decay check
    dcy_pt = mp.Vector3(support_thick/2)
    dcy_cn = {'s':mp.Ez, 'p':mp.Ey}[pol]

    #Run sim until decays
    sim.run(until_after_sources=mp.stop_when_fields_decayed( \
        dt=decay_dt, pt=dcy_pt, c=dcy_cn, decay_by=decay_by))

    #Synchronize magnetic fields
    sim.fields.synchronize_magnetic_fields()

    #Output DFT fields
    sim.output_dft(dft_obj, dft_name)

    #Get metadata
    if get_meta:

        #Get coordinates
        x,y,z,w = sim.get_array_metadata(dft_cell=dft_obj)

        #Save metadata
        if mp.am_master():

            #Prefix
            pre = f"{save_dir}/{['', 'vac-'][int(is_vac)]}"

            #Save coordinates
            with h5py.File(f'{pre}meta.h5', 'w') as f:
                f.create_dataset('xx', data=x)
                f.create_dataset('yy', data=y)
                f.create_dataset('waves', data=waves)
                f.create_dataset('gap_width', data=gap_width)
                f.create_dataset('device_bottom', data=device_bottom)
                f.create_dataset('support_bottom', data=support_bottom)

        #Cleanup
        del x, y, z, w

    #Reset meep
    sim.reset_meep()

############################################
####	Loop over angles and run sim ####
############################################

#Get meta once
get_meta = True

for pol in ['s','p']:
    for is_vac in [False, True]:

        #Run sim
        run_sim(pol, is_vac, get_meta)

        #Clear flag
        get_meta = False
