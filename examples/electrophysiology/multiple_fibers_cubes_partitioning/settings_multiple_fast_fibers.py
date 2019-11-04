# Multiple 1D fibers (monodomain), biceps geometry
# This is similar to the fibers_emg example, but without EMG.
# To see all available arguments, execute: ./multiple_fibers settings_multiple_fibers_cubes_partitioning.py -help
#
# if fiber_file=cuboid.bin, it uses a small cuboid test example
#
# You have to set n_subdomains such that it matches the number of processes, e.g. 2x2x1 = 4 processes.
# Decomposition is in x,y,z direction, the fibers are aligned with the z axis.
# E.g. --n_subdomains 2 2 1 which is 2x2x1 means no subdivision per fiber, 
# --n_subdomains 8 8 4 means every fiber will be subdivided to 4 processes and all fibers will be computed by 8x8 processes.
#
# Example with 4 processes and end time 5, and otherwise default parameters:
#   mpirun -n 4 ./multiple_fibers settings_multiple_fibers_cubes_partitioning.py --n_subdomains 2 2 1 --end_time=5.0
#
# Three files contribute to the settings:
# A lot of variables are set by the helper.py script, the variables and their values are defined in variables.py and this file
# creates the composite config that is needed by opendihu.
# You can provided parameter values in a custom_variables.py file in the variables subfolder. (Instead of custom_variables.py you can choose any filename.)
# This custom variables file should be the next argument on the command line after settings_fibers_emg.py, e.g.:
#
#  ./multiple_fibers settings_multiple_fibers_cubes_partitioning.py custom_variables.py --n_subdomains 1 1 1 --end_time=5.0

import sys, os
import timeit
import argparse
import importlib

# parse rank arguments
rank_no = (int)(sys.argv[-2])
n_ranks = (int)(sys.argv[-1])

# add variables subfolder to python path where the variables script is located
script_path = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, script_path)
sys.path.insert(0, os.path.join(script_path,'variables'))

import variables              # file variables.py, defined default values for all parameters, you can set the parameters there  
from create_partitioned_meshes_for_settings import *   # file create_partitioned_meshes_for_settings with helper functions about own subdomain
import settings_multiple_fibers_cubes_partitioning

config = settings_multiple_fibers_cubes_partitioning.config
