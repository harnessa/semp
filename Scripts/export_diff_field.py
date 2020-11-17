"""
export_diff_field.py

Author: Anthony Harness
Affiliation: Princeton University
Created on: 11-16-2020
Package: SEMP

Description: Export difference field from SEMP simulation for use in greypixel
License: Refer to $pkg_home_dir/LICENSE
"""

import semp
import numpy as np

params = {
    ### Saving ###
    # 'do_save':          True,
    # 'save_dir':         ,

    ### Sim ###
    # 'session':          'corner_test',
    'session':          'test',
    # 'ext':              'not_full',
    'ext':              'full',
    'polarization':     's',
}

alz = semp.analysis.Analyzer(params)
alz.export_difference_field()
