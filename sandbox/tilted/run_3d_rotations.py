import numpy as np
import meep as mp
import meep.materials as mat_lib
import h5py
import os

############################################
####	User Input ####
############################################

#Numerics
resolution = 40
dpml = 4
dpad = 4
seam_dark = 15
seam_lite = 45

#Wafer
wafer_thick = 2.25
skin_thick = 0.3
wafer_material = 'cSi'
skin_material = 'Al'

#Scallops
scallop_list = [
    [(0.205, 0.324), (0.762, 0.810)],
    [(0.148, 1.157), (0.762, 0.824)],
    [(0.019, 1.900), (0.762, 0.767)]
]

#Angles to run
inn_ang_max = 10
out_ang_max = 10
nangs = 11

#Directory
# base_dir = '/home/aharness/Research/Optics_Modeling/Semp_Results/test'
# base_dir = '/home/aharness/Research/Optics_Modeling/Semp_Results/tilted_runs'
base_dir = '/scratch/network/aharness/Semp_Results/tilted_runs/scallops_3D_r40'

############################################
####	Other params ####
############################################

#Decay
decay_dt = 50         # Time after decay to
decay_by = 1e-6       # Decay amount
waves = np.array([0.641, 0.660, 0.699, 0.725])

#Force to be odd (so that we get 0 angle)
nangs += (nangs + 1) % 2

#Make sure x big enough to get full seam
xseam_pad = max(2*dpad, np.ceil((seam_lite + 0.5*seam_dark)*np.sin(np.radians(inn_ang_max))))

#Make sure seam_lite is big enough for largest angle
seam_lite = max(seam_lite, seam_dark + wafer_thick/np.sin(np.radians(inn_ang_max)))

#Derived
freqs = 1./waves
fcen0 = freqs.mean()
dfreq = max(2.*freqs.ptp(), 0.35)
inn_angs = np.linspace(-inn_ang_max, inn_ang_max, nangs)
out_angs = np.linspace(0, out_ang_max, nangs)

mat_lib.metal = mp.metal
waf_mat = getattr(mat_lib, wafer_material)
skn_mat = getattr(mat_lib, skin_material)

lx_np = wafer_thick + 2*xseam_pad
ly_np = dpad*2 + seam_lite + seam_dark
cell_size = mp.Vector3(dpml*2 + lx_np, dpml*2 + ly_np)

############################################
####	Simulation Function ####
############################################

def run_sim(inn_ang, out_ang, pol, is_vac, get_meta):

    #Print
    if mp.am_master():
        print(f'\nRunning: {inn_ang}, {out_ang}, {pol}, {is_vac}\n')

    #Run specific
    inn_rot_angle = np.radians(inn_ang)
    out_rot_angle = np.radians(out_ang)
    src_comp = getattr(mp, {'s': 'Ez', 'p': 'Ey'}[pol])

    #Save dir
    sgn = ["n", "p"][int((np.sign(inn_rot_angle)+1)/2)]
    save_dir = f'{base_dir}/{sgn}{abs(inn_ang)*100:.0f}_{abs(out_ang)*100:.0f}'
    if mp.am_master():
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)

    #Source
    sim_src = mp.GaussianSource(fcen0, fwidth=dfreq, is_integrated=True)
    sources = [mp.Source(sim_src, component=src_comp,
                center=mp.Vector3(x=-0.5*lx_np), \
                size=mp.Vector3(y=cell_size.y))]

    #PML
    if pol == 's' or is_vac:
        pml_layers = [mp.PML(thickness=dpml)]
    else:
        #Don't put PML on wafer side -- blows up
        pml_layers = [mp.PML(thickness=dpml, direction=mp.X), \
            mp.PML(thickness=dpml, direction=mp.Y, side=-1), \
            mp.Absorber(thickness=dpml, direction=mp.Y, side=1)]

    #Geometry
    geometry = []
    if not is_vac:

        #Rotate function
        rot_func = lambda vec: vec.rotate(mp.Vector3(z=1), inn_rot_angle)

        #Wafer
        waf_sy = dpad+dpml+seam_dark
        waf_cc = mp.Vector3(y=-0.5*cell_size.y + 0.4*waf_sy)
        e1 = rot_func(mp.Vector3(x=1))
        e2 = rot_func(mp.Vector3(y=1))
        wafer = mp.Block(material=waf_mat, size=mp.Vector3(wafer_thick, waf_sy, mp.inf),
            center=mp.Vector3(y=waf_cc.y), e1=e1, e2=e2)

        #Skin (Rotate about center of wafer)
        skn_cc = mp.Vector3(-wafer_thick/2 - skin_thick/2, waf_cc.y)
        skn_cc = rot_func(skn_cc - waf_cc) + waf_cc
        skin = mp.Block(material=skn_mat ,size=mp.Vector3(skin_thick, waf_sy, mp.inf),
            center=mp.Vector3(skn_cc.x, skn_cc.y), e1=e1, e2=e2)

        #Get edge of wafer
        edge_y = wafer.center.y + wafer.size.y/2
        #Add to list
        geometry += [wafer, skin]

        #Scallops
        if len(scallop_list) > 0:
            #Get top point of wafer
            waf_p0 = mp.Vector3(-wafer_thick/2, edge_y).rotate(mp.Vector3(z=1), inn_rot_angle)

            #Loop through scallop list
            for scallop in scallop_list:
                #Get current size and center
                cur_cen, cur_sze = scallop
                scl_sze = mp.Vector3(cur_sze[1], cur_sze[0], mp.inf)
                #Rotate center about center of wafer
                scl_cc = mp.Vector3(-wafer_thick/2 + cur_cen[1], edge_y + cur_cen[0])
                scl_cc = rot_func(scl_cc - waf_cc) + waf_cc
                #Add air ellipsoid
                geometry += [mp.Ellipsoid(material=mp.air, size=scl_sze, \
                    center=scl_cc, e1=e1, e2=e2)]

    #Kpoint
    kpoint = mp.Vector3(z=np.sin(out_rot_angle)).scale(fcen0)

    #Build simulation
    sim = mp.Simulation(force_complex_fields=True,
        resolution=resolution, ensure_periodicity=False,
        cell_size=cell_size, boundary_layers=pml_layers, sources=sources,
        geometry=geometry, k_point=kpoint)

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
                f.create_dataset('inn_rot_angle', data=inn_rot_angle)
                f.create_dataset('out_rot_angle', data=out_rot_angle)
                f.create_dataset('wafer_center', data=[wafer.center.x, wafer.center.y])
                f.create_dataset('wafer_thick', data=wafer_thick)
                f.create_dataset('skin_thick', data=skin_thick)
                f.create_dataset('resolution', data=resolution)
                f.create_dataset('dpml_dpad', data=[dpml, dpad])
                f.create_dataset('seams__dark_lite', data=[seam_dark, seam_lite])
                f.create_dataset('edge_y', data=edge_y)

        #Cleanup
        del x, y, z, w

    #Reset meep
    sim.reset_meep()

############################################
####	Loop over angles and run sim ####
############################################

for inn_ang in inn_angs:
    for out_ang in out_angs:
        #Get meta for each angle
        get_meta = True

        for pol in ['s','p']:
            for is_vac in [False, True]:

                #Run sim
                run_sim(inn_ang, out_ang, pol, is_vac, get_meta)

                #Clear flag
                get_meta = False
