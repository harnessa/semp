import semp
import matplotlib.pyplot as plt;plt.ion()

for pol in ['s', 'p']:

    MEEP_params = {
        ### Lab Properties  ###
        'polarization':     pol,
        'wave':             0.641,

        ### Mask Properties ###
        'sim_geometry':     'vacuum',
        'is_sommerfeld':    True,

        ### Numerics ###
        'resolution':       20,
        'pml_all':          4,
        'pad_all':          4,
        'n_periods':        80,

        'obs_distance': 0,

    }

    PROP_params = {'do_save': True, 'is_movie': True}

    #Run simulation
    prop = semp.Propagator(MEEP_params, PROP_params)
    prop.run_sim()
