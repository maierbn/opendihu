# Linear elasticity
# Use paraview and the Warp Filter for visualization.
# Add a Glyph filter with arrows for field "-rhsNeumannBC", this is the applied traction.

import numpy as np
import sys, os
import argparse
import distutils.util
sys.path.insert(0, '..')
sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
import variables              # file variables.py, defined default values for all parameters, you can set the parameters there
from create_partitioned_meshes_for_settings import *   # file create_partitioned_meshes_for_settings with helper functions about own subdomain

# define command line arguments
mbool = lambda x:bool(distutils.util.strtobool(x))   # function to parse bool arguments
parser = argparse.ArgumentParser(description='3d_hyperelasticity')
parser.add_argument('--yforce',                help='Applied force in y direction.',      type=float,   default=-0.4)
parser.add_argument('--zforce',                help='Applied force in z direction.',      type=float,   default=-3)

# parse command line arguments and assign values to variables module
args, other_args = parser.parse_known_args(args=sys.argv[:-2], namespace=variables)
if len(other_args) != 0:
    print("Warning: These arguments were not parsed by the settings python file\n  " + "\n  ".join(other_args), file=sys.stderr)

# muscle mesh
# -----------
# stride for sampling the 3D elements from the fiber data
# here any number is possible
sampling_stride_x = 7
sampling_stride_y = 7
sampling_stride_z = 500

sampling_stride_x = 1
sampling_stride_y = 1
sampling_stride_z = 95

# input mesh file
input_directory = os.path.join(os.environ["OPENDIHU_HOME"], "examples/electrophysiology/input")
fiber_file = input_directory+"/left_biceps_brachii_13x13fibers.bin"

load_fiber_data = True             # If the fiber geometry data should be loaded completely in the python script. If True, this reads the binary file and assigns the node positions in the config. If False, the C++ code will read the binary file and only extract the local node positions. This is more performant for highly parallel runs.

# create the partitioning using the script in create_partitioned_meshes_for_settings.py
result = create_partitioned_meshes_for_settings(
    variables.n_subdomains_x, variables.n_subdomains_y, variables.n_subdomains_z, 
    fiber_file, load_fiber_data,
    sampling_stride_x, sampling_stride_y, sampling_stride_z, generate_linear_3d_mesh=True, generate_quadratic_3d_mesh=True)

#parse result
[variables.meshes, variables.own_subdomain_coordinate_x, variables.own_subdomain_coordinate_y, variables.own_subdomain_coordinate_z, variables.n_fibers_x, variables.n_fibers_y, variables.n_points_whole_fiber] = result

# debugging output
node_positions = variables.meshes["3Dmesh_quadratic"]["nodePositions"]

[mx, my, mz] = variables.meshes["3Dmesh_quadratic"]["nPointsGlobal"]
[nx, ny, nz] = variables.meshes["3Dmesh_quadratic"]["nElements"]
print("node_positions: {}, n elements: {},{},{}, n nodes: {},{},{}".format(len(node_positions),nx,ny,nz,mx,my,mz))
print("applied force: y: {}, z: {}".format(variables.yforce, variables.zforce))
#dirichlet_bc = {}
#neumann_bc = []

# set Dirichlet BC at top nodes for linear elasticity problem, fix muscle at top
elasticity_dirichlet_bc = {}
for j in range(my):
  for i in range(mx):
    elasticity_dirichlet_bc[(mz-1)*mx*my + j*mx + i] = [None,None,0.0]
  
# fix edge
for i in range(mx):
  elasticity_dirichlet_bc[(mz-1)*mx*my + 0*mx + i] = [None,0.0,0.0]
  
# fix corner completely
elasticity_dirichlet_bc[(mz-1)*mx*my + 0] = [0.0,0.0,0.0]

# traction force, given by command line arguments --yforce and --zforce
bottom_traction = [0.0,variables.yforce,variables.zforce]    # [N]


# Neumann BC at bottom nodes, traction downwards
elasticity_neumann_bc = [{"element": 0*nx*ny + j*nx + i, "constantVector": bottom_traction, "face": "2-"} for j in range(ny) for i in range(nx)]

# material parameters
# --------------------
bulk_modulus = 3.9  # [10 kPa] 39 kPa   # https://www.researchgate.net/publication/230248067_Bulk_Modulus
shear_modulus = 4.8 # [10 kPa] 48 kPa   # https://onlinelibrary.wiley.com/doi/full/10.1002/mus.24104

# function to postprocess the output
def postprocess(result):
  result = result[0]
  #print(result)
  
  current_time = result["currentTime"]
  timestep_no = result["timeStepNo"]
  
  field_variables = result["data"]
  # field_variables[0] is the geometry
  # field_variables[1] is the displacements u
  # field_variables[2] is the velocities v
  # field_variables[3] is the PK2-Stress (Voigt)
  u_components = field_variables[1]["components"]
  # stress_components contains the symmetric 3x3 material stress tensor, in this order: S_11, S_22, S_33, S_12, S_13, S_23

  ux = u_components[0]["values"]
  uy = u_components[1]["values"]
  uz = u_components[2]["values"]
  
  print("displacements: {},{},{}".format(ux[0],uy[0],uz[0]))
  with open("displacements.csv", "a") as f:
    f.write("linear;{};{}\n".format(variables.zforce,uz[0]))

config = {
  "scenarioName":                   "linear_elasticity_3d",          # scenario name for the log file
  "logFormat":                      "csv",                           # "csv" or "json", format of the lines in the log file, csv gives smaller files
  "solverStructureDiagramFile":     "out/solver_structure.txt",          # output file of a diagram that shows data connection between solvers
  "mappingsBetweenMeshesLogFile":   "out/mappings_between_meshes.txt",   # log file for mappings between meshes
  
  "Meshes": variables.meshes,
  "FiniteElementMethod": {
    
    # mesh
    "meshName": "3Dmesh_quadratic",                 # muscle mesh with quadratic Lagrange ansatz functions
    "inputMeshIsGlobal": True,           # boundary conditions are specified in global numberings
    
    "outputInterval":     1.0,                      # the timestep interval in which to print the current time step to output
    "prefactor": 1,                                 # prefactor c of the equation c*Δu = f
    "dirichletBoundaryConditions":  elasticity_dirichlet_bc,   # dirichlet boundary conditions
    "dirichletOutputFilename":      None,           # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
    "neumannBoundaryConditions":    elasticity_neumann_bc,     # neumann boundary conditions
    "divideNeumannBoundaryConditionValuesByTotalArea": True,  # if the neumann boundary condition vectors should be divided by the total surface area where surface loads are applied, this allows to specify the total force that acts on the surface. If set to False (default), the given traction is a per-surface quantity.
    "slotName":           "",                       # slot name of the solution variable
    
    # solver
    "relativeTolerance":  1e-15,                    # relative tolerance for the solver, this relates to the current residual norm divided by the norm of the rhs
    "absoluteTolerance":  1e-10,                    # absolute tolerance of the residual        
    "solverType":         "gmres",                  # which PETSc linear system solver to use
    "preconditionerType": "none",                 # which PETSc preconditioner to use
    "maxIterations":      1e4,                      # maximum number of iterations of the linear solver
    "dumpFormat":         "matlab",                 # "default", "ascii", "matlab", output format type for system matrix and right hand side
    "dumpFilename":       "",                       # filename for output of system matrix and right hand side after every solve, set to "" to disable
    
    # material parameters
    "bulkModulus":        bulk_modulus,                      # bulk modulus, K, material parameter for compressibility
    "shearModulus":       shear_modulus,                      # shear modulus, μ
    
    "OutputWriter": [
      {"format": "Paraview", "outputInterval": 1, "filename": "out/3d", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
      {"format": "PythonCallback", "callback": postprocess, "outputInterval": 1, "binary":True, "onlyNodalValues":True, "fileNumbering": "incremental"},
    ]
  },
}
