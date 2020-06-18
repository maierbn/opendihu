# biceps for chaste solver
#

import numpy as np
import pickle
import sys
sys.path.insert(0, '..')
import variables              # file variables.py, defined default values for all parameters, you can set the parameters there
from create_partitioned_meshes_for_settings import *   # file create_partitioned_meshes_for_settings with helper functions about own subdomain


# input mesh file
fiber_file = "../../../electrophysiology/input/13x13fibers.bin"

load_fiber_data = True             # If the fiber geometry data should be loaded completely in the python script. If True, this reads the binary file and assigns the node positions in the config. If False, the C++ code will read the binary file and only extract the local node positions. This is more performant for highly parallel runs.

# partitioning
# ------------
# this has to match the total number of processes
n_subdomains_x = 1
n_subdomains_y = 1
n_subdomains_z = 1

# stride for sampling the 3D elements from the fiber data
# here any number is possible
sampling_stride_x = 7
sampling_stride_y = 7
sampling_stride_z = 1000

#sampling_stride_x = 1
#sampling_stride_y = 1
#sampling_stride_z = 100

# create the partitioning using the script in create_partitioned_meshes_for_settings.py
result = create_partitioned_meshes_for_settings(
    n_subdomains_x, n_subdomains_y, n_subdomains_z, 
    fiber_file, load_fiber_data,
    sampling_stride_x, sampling_stride_y, sampling_stride_z, True)

#parse result
[variables.meshes, variables.own_subdomain_coordinate_x, variables.own_subdomain_coordinate_y, variables.own_subdomain_coordinate_z, variables.n_fibers_x, variables.n_fibers_y, variables.n_points_whole_fiber] = result

node_positions = variables.meshes["3Dmesh_quadratic"]["nodePositions"]
#node_positions = variables.meshes["3Dmesh"]["nodePositions"]

"""
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


x_list = [x for [x,y,z] in node_positions]
y_list = [y for [x,y,z] in node_positions]
z_list = [z for [x,y,z] in node_positions]
  
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x_list,y_list,z_list,c='b',s=30)

plt.show()
"""

config = {
  "scenarioName": "3d_muscle_chaste",
  "logFormat":    "csv",     # "csv" or "json", format of the lines in the log file, csv gives smaller files
  "solverStructureDiagramFile":     "solver_structure.txt",     # output file of a diagram that shows data connection between solvers
  "Meshes": variables.meshes,
  "QuasiStaticNonlinearElasticitySolverChaste": {
    "maximumActiveStress": 1.0,       # dummy value
    "strainScalingCurveWidth": 1.0,    # dummy value
    "meshName": "3Dmesh_quadratic",
    "OutputWriters": [
      {"format": "PythonFile", "outputInterval": 1, "filename": "out/mesh", "binary": True}
    ]
  }
}
