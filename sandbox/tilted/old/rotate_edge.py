import numpy as np
import meep as mp
import h5py

import matplotlib.pyplot as plt;plt.ion()

############################################
####	User Input ####
############################################

is_vac = [False, True][1]

resolution = 10

rot_angle = np.radians(0)

angs = [0, 10, -10][:1]

for ang in angs:

    rot_angle = np.radians(ang)

    pol = ['s','p'][0]

    dpml = 2
    dpad = 3
    gap = 5
    thick = 1

    wave = 0.6
    df = 0.02  # turn-on bandwidth

    ############################################
    ############################################

    #Derived
    fcen = 1/wave
    lx_np = dpad*2 + thick
    ly_np = dpad*2 + gap
    cell_size = mp.Vector3(dpml*2 + lx_np, dpml*2 + ly_np)
    src_comp = getattr(mp, {'s': 'Ez', 'p': 'Ey'}[pol])

    #Kpoint
    theta = np.radians(5)
    kpoint = mp.Vector3(z=np.sin(theta)).scale(fcen)

    def pw_amp(k, x0):
        def _pw_amp(x):
            return np.exp(1j * 2*np.pi*k.dot(x + x0))
        return _pw_amp

    #Source
    sim_src = mp.GaussianSource(fcen, fwidth=df, is_integrated=True)
    # sim_src = mp.ContinuousSource(fcen, fwidth=df, is_integrated=True)
    src_pt = mp.Vector3(x=-0.5*lx_np)
    sources = [mp.Source(sim_src, component=src_comp,
                center=src_pt, \
                size=mp.Vector3(y=cell_size.y),
                amp_func=pw_amp(kpoint, src_pt))]

    #PML
    if pol == 's' or is_vac:
        pml_layers = [mp.PML(thickness=dpml)]
    else:
        pml_layers = [mp.PML(thickness=dpml,  direction=mp.X), \
            mp.PML(thickness=dpml, direction=mp.Y, side=-1),
            mp.Absorber(thickness=dpml, direction=mp.Y, side=1)]

    #Geometry
    geometry = []
    if not is_vac:
        waf_sy = dpad+dpml
        e1 = mp.Vector3(x=1).rotate(mp.Vector3(z=1), rot_angle)
        e2 = mp.Vector3(y=1).rotate(mp.Vector3(z=1), rot_angle)
        wafer = mp.Block(material=mp.metal, size=mp.Vector3(thick, waf_sy, mp.inf),
            center=mp.Vector3(y=-0.5*cell_size.y + 0.4*waf_sy), e1=e1, e2=e2)
        geometry += [wafer]

    #Build simulation
    sim = mp.Simulation(force_complex_fields=True,
        resolution=resolution, ensure_periodicity=False,
        cell_size=cell_size, boundary_layers=pml_layers, sources=sources,
        geometry=geometry, k_point=kpoint)

    ############################################
    ############################################

    #Add DFT fields
    vol = mp.Volume(size=mp.Vector3(lx_np, ly_np))
    all_flds = [mp.Ex, mp.Ey, mp.Ez, mp.Hx, mp.Hy, mp.Hz]
    dft_obj = sim.add_dft_fields(all_flds, [fcen], where=vol)

    #Run
    sim.run(until_after_sources=mp.stop_when_fields_decayed(50, src_comp, mp.Vector3(-0.5*cell_size.x+dpml), 1e-6))
    # sim.run(until=400)
    # sim.plot2D(fields=src_comp, output_plane=vol)

    # fld = sim.get_array(src_comp, vol=vol)
    fld = sim.get_dft_array(dft_obj, src_comp, 0)
    x,y,z,w = sim.get_array_metadata(dft_cell=dft_obj)

    # if mp.am_master():
    #     with h5py.File(f'./{pol}_{["n", "p"][int((np.sign(rot_angle)+1)/2)]}{abs(np.degrees(rot_angle)):.0f}.h5', 'w') as f:
    #         f.create_dataset('field', data=fld)
    #         f.create_dataset('x', data=x)
    #         f.create_dataset('y', data=y)
    #         f.create_dataset('rot_angle', data=rot_angle)

    ff = ['Ex', 'Ey', 'Ez', 'Hx', 'Hy', 'Hz']

    for f in ff:
        fld = sim.get_dft_array(dft_obj, getattr(mp, f), 0)
        if fld.size < 2:
            continue
        print(f, abs(fld).max())
        # plt.figure()
        # plt.imshow(abs(fld))
        # plt.title(f)
        #

    # plt.figure()
    # plt.imshow(np.abs(fld))
    #
    # plt.figure()
    # plt.imshow(np.angle(fld))
    breakpoint()
    # # fig, axes = plt.subplots(1,2)
    # # axes[0].imshow(np.abs(fld).T, origin='lower')
    # # axes[1].imshow(np.angle(fld).T, origin='lower')
    #
    # # #Get field
    # # fld = sim.get_dft_array(dft_obj, src_comp, 0)
    #
    #
    # breakpoint()
