"""
__init__.py

Author: Anthony Harness
Affiliation: Princeton University
Created on: 10-06-2020
Package: SEMP

Description: __init__ package for the SEMP python package
License: Refer to $pkg_home_dir/LICENSE
"""

import numpy as np
import os

#####################
#####   MPI #####
#####################

try:
    from mpi4py import MPI
    mpi_rank = MPI.COMM_WORLD.rank      # processor ID number, from 0 up to size
    mpi_size = MPI.COMM_WORLD.size      # total number of processors running
    mpi_barrier = MPI.COMM_WORLD.Barrier
    has_mpi = True
except ImportError:
    mpi_rank = 0
    mpi_size = 1
    mpi_barrier = lambda : None
    has_mpi = False
zero_rank = mpi_rank == 0

#####################
#####   Directories #####
#####################

pkg_home_dir = os.getenv("SEMP")

if pkg_home_dir is None:
    if zero_rank:
        print("\n*** Cannot Find Environment Variable pointing to SEMP home! ***\n")
        print("* Please set environment variable $SEMP pointing to directory where semp/setup.py is located *")
    import sys
    sys.exit()

results_dir = "%s/Research/Optics_Modeling/Semp_Results"%os.getenv("HOME")
ext_data_dir = f"{pkg_home_dir}/External_Data"
int_data_dir = f"{pkg_home_dir}/Internal_Data"

#####################
#####   Modules #####
#####################

import semp.utils
from .propagator import Propagator
import semp.simulation
import semp.analysis
