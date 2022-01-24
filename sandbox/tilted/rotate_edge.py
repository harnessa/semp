import numpy as np
import meep as mp
import h5py

import matplotlib.pyplot as plt;plt.ion()

############################################
####	User Input ####
############################################

is_vac = [False, True][0]

resolution = 10

rot_angle = np.radians(0)

angs = [0, 10, -10][1:]

for ang in angs:

    rot_angle = np.radians(ang)

    pol = ['s','p'][0]

    dpml = 2
    dpad = 3
    gap = 4
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

    #Source
    sim_src = mp.GaussianSource(fcen, fwidth=df, is_integrated=True)
    # sim_src = mp.ContinuousSource(fcen, fwidth=df, is_integrated=True)
    sources = [mp.Source(sim_src, component=src_comp,
                center=mp.Vector3(x=-0.5*lx_np), \
                size=mp.Vector3(y=cell_size.y))]

    #PML
    pml_layers = [mp.PML(thickness=dpml)]

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
        geometry=geometry, k_point=mp.Vector3())

    ############################################
    ############################################

    #Add DFT fields
    vol = mp.Volume(size=mp.Vector3(lx_np, ly_np))
    dft_obj = sim.add_dft_fields([src_comp], [fcen], where=vol)

    vol2 = mp.Volume(center=mp.Vector3(x=thick/2), size=mp.Vector3(1, ly_np))

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

    plt.figure()
    plt.imshow(abs(fld))
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
