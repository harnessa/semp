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
seam_dark = 15
seam_lite = 35

#Wafer
wafer_thick = 2.25
skin_thick = 0.3
wafer_material = 'cSi'
skin_material = 'Al'

#Angles to run
ang_max = 5
nangs = 10

############################################
####	Other params ####
############################################

#Decay
decay_dt = 50         # Time after decay to
decay_by = 1e-3       # Decay amount
waves = np.array([0.641, 0.660, 0.699, 0.725])
base_dir = '/home/aharness/Research/Optics_Modeling/Semp_Results/tilted_runs'

#Derived
freqs = 1./waves
fcen0 = freqs.mean()
dfreq = max(2.*freqs.ptp(), 0.35)
angs = np.linspace(-ang_max, ang_max, nangs)

mat_lib.metal = mp.metal
waf_mat = getattr(mat_lib, wafer_material)
skn_mat = getattr(mat_lib, skin_material)

lx_np = dpad*2 + wafer_thick
ly_np = dpad*2 + seam_lite + seam_dark
cell_size = mp.Vector3(dpml*2 + lx_np, dpml*2 + ly_np)

############################################
####	Simulation Function ####
############################################

def run_sim(ang, pol, is_vac, get_meta):

    #Run specific
    rot_angle = np.radians(ang)
    src_comp = getattr(mp, {'s': 'Ez', 'p': 'Ey'}[pol])

    #Save dir
    sgn = ["n", "p"][int((np.sign(rot_angle)+1)/2)]
    save_dir = f'{base_dir}/{sgn}{abs(ang)*100:.0f}'
    if mp.am_master():
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)

    #Source
    sim_src = mp.GaussianSource(fcen0, fwidth=dfreq, is_integrated=True)
    sources = [mp.Source(sim_src, component=src_comp,
                center=mp.Vector3(x=-0.5*lx_np), \
                size=mp.Vector3(y=cell_size.y))]

    #PML
    pml_layers = [mp.PML(thickness=dpml)]

    #Geometry
    geometry = []
    if not is_vac:
        waf_sy = dpad+dpml+seam_dark
        waf_cy = -0.5*cell_size.y + 0.4*waf_sy
        e1 = mp.Vector3(x=1).rotate(mp.Vector3(z=1), rot_angle)
        e2 = mp.Vector3(y=1).rotate(mp.Vector3(z=1), rot_angle)
        wafer = mp.Block(material=waf_mat, size=mp.Vector3(wafer_thick, waf_sy, mp.inf),
            center=mp.Vector3(y=waf_cy), e1=e1, e2=e2)
        skin = mp.Block(material=skn_mat ,size=mp.Vector3(skin_thick, waf_sy, mp.inf),
            center=mp.Vector3(-wafer_thick/2 - skin_thick/2, waf_cy), e1=e1, e2=e2)
        geometry += [wafer, skin]

    #Build simulation
    sim = mp.Simulation(force_complex_fields=True,
        resolution=resolution, ensure_periodicity=False,
        cell_size=cell_size, boundary_layers=pml_layers, sources=sources,
        geometry=geometry, k_point=mp.Vector3())

    ############################################
    ############################################

    #Set output directory + prefix
    sim.use_output_directory(save_dir)
    sim.filename_prefix = ['', 'vac'][int(is_vac)]
    dft_name = f'{save_dir}/{["", "vac-"][int(is_vac)]}fields_{pol}'

    #Get fields to output
    fld_names = {'s':[mp.Ez, mp.Hy], 'p':[mp.Hz, mp.Ey]}[pol]
    far_names = {'s':['Ez', 'Hy'], 'p':['Hz', 'Ey']}[pol]

    #Add DFT fields
    vol = mp.Volume(size=mp.Vector3(lx_np, ly_np))
    dft_obj = sim.add_dft_fields(fld_names, freqs, where=vol)

    #Decay check
    dcy_pt = mp.Vector3(wafer_thick/2, ly_np/2)
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
                f.create_dataset('rot_angle', data=rot_angle)
                f.create_dataset('wafer_center', data=[wafer.center.x, wafer.center.y])

        #Cleanup
        del x, y, z, w

    #Reset meep
    sim.reset_meep()

############################################
####	Loop over angles and run sim ####
############################################

for ang in angs:
    #Get meta for each angle
    get_meta = True

    for pol in ['s','p']:
        for is_vac in [False, True]:

            #Run sim
            run_sim(ang, pol, is_vac, get_meta)

            #Clear flag
            get_meta = False
