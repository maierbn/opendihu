# biceps
#

import numpy as np
import pickle
import argparse
import distutils.util
import sys, os
sys.path.insert(0, '..')
sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
import variables              # file variables.py, defined default values for all parameters, you can set the parameters there
from create_partitioned_meshes_for_settings import *   # file create_partitioned_meshes_for_settings with helper functions about own subdomain

own_rank_no = (int)(sys.argv[-2])
n_ranks = (int)(sys.argv[-1])

# define command line arguments
mbool = lambda x:bool(distutils.util.strtobool(x))   # function to parse bool arguments
parser = argparse.ArgumentParser(description='3d_hyperelasticity')
parser.add_argument('--scenario_name',         help='The name to identify this run in the log.',      default="3d_muscle")
parser.add_argument('--njacobi',               help='Frequency of Jacobi recomputation.', type=int,   default=5)
parser.add_argument('--yforce',                help='Applied force in y direction.',      type=float,   default=-0.4)
parser.add_argument('--zforce',                help='Applied force in z direction.',      type=float,   default=-3)
parser.add_argument('--n_subdomains', nargs=3, help='Number of subdomains in x,y,z direction.',             type=int)

# parse command line arguments and assign values to variables module
args, other_args = parser.parse_known_args(args=sys.argv[:-2], namespace=variables)
if len(other_args) != 0 and own_rank_no == 0:
    print("Warning: These arguments were not parsed by the settings python file\n  " + "\n  ".join(other_args), file=sys.stderr)


# input mesh file
input_directory = os.path.join(os.environ["OPENDIHU_HOME"], "examples/electrophysiology/input")
fiber_file = input_directory+"/left_biceps_brachii_13x13fibers.bin"

load_fiber_data = True             # If the fiber geometry data should be loaded completely in the python script. If True, this reads the binary file and assigns the node positions in the config. If False, the C++ code will read the binary file and only extract the local node positions. This is more performant for highly parallel runs.


# initialize some dependend variables
if variables.n_subdomains is not None:
  variables.n_subdomains_x = variables.n_subdomains[0]
  variables.n_subdomains_y = variables.n_subdomains[1]
  variables.n_subdomains_z = variables.n_subdomains[2]
  
variables.n_subdomains = variables.n_subdomains_x*variables.n_subdomains_y*variables.n_subdomains_z

# partitioning
# ------------

# automatically initialize partitioning if it has not been set
if n_ranks != variables.n_subdomains:
  
  # create all possible partitionings to the given number of ranks
  optimal_value = n_ranks**(1/3)
  possible_partitionings = []
  for i in range(1,n_ranks+1):
    for j in range(1,n_ranks+1):
      if i*j <= n_ranks and n_ranks % (i*j) == 0:
        k = (int)(n_ranks / (i*j))
        performance = (k-optimal_value)**2 + (j-optimal_value)**2 + 1.1*(i-optimal_value)**2
        possible_partitionings.append([i,j,k,performance])
        
  # if no possible partitioning was found
  if len(possible_partitionings) == 0:
    if rank_no == 0:
      print("\n\n\033[0;31mError! Number of ranks {} does not match given partitioning {} x {} x {} = {} and no automatic partitioning could be done.\n\n\033[0m".format(n_ranks, variables.n_subdomains_x, variables.n_subdomains_y, variables.n_subdomains_z, variables.n_subdomains_x*variables.n_subdomains_y*variables.n_subdomains_z))
    quit()
    
  # select the partitioning with the lowest value of performance which is the best
  lowest_performance = possible_partitionings[0][3]+1
  for i in range(len(possible_partitionings)):
    if possible_partitionings[i][3] < lowest_performance:
      lowest_performance = possible_partitionings[i][3]
      variables.n_subdomains_x = possible_partitionings[i][0]
      variables.n_subdomains_y = possible_partitionings[i][1]
      variables.n_subdomains_z = possible_partitionings[i][2]

# stride for sampling the 3D elements from the fiber data
# here any number is possible
sampling_stride_x = 7
sampling_stride_y = 7
sampling_stride_z = 500

sampling_stride_x = 1
sampling_stride_y = 1
sampling_stride_z = 95

# create the partitioning using the script in create_partitioned_meshes_for_settings.py
result = create_partitioned_meshes_for_settings(
    variables.n_subdomains_x, variables.n_subdomains_y, variables.n_subdomains_z, 
    fiber_file, load_fiber_data,
    sampling_stride_x, sampling_stride_y, sampling_stride_z, generate_linear_3d_mesh=True, generate_quadratic_3d_mesh=True)

#parse result
[variables.meshes, variables.own_subdomain_coordinate_x, variables.own_subdomain_coordinate_y, variables.own_subdomain_coordinate_z, variables.n_fibers_x, variables.n_fibers_y, variables.n_points_whole_fiber] = result

node_positions = variables.meshes["3Dmesh_quadratic"]["nodePositions"]
#node_positions = variables.meshes["3Dmesh"]["nodePositions"]

print("applied force: y: {}, z: {}".format(variables.yforce, variables.zforce))

# material parameters
# --------------------
# quantities in mechanics unit system
rho = 10                    # [1e-4 kg/cm^3] density of the muscle (density of water)

# Mooney-Rivlin parameters [c1,c2,b,d] of c1*(Ibar1 - 3) + c2*(Ibar2 - 3) + b/d (λ - 1) - b*ln(λ)
# Heidlauf13: [6.352e-10 kPa, 3.627 kPa, 2.756e-5 kPa, 43.373] = [6.352e-11 N/cm^2, 3.627e-1 N/cm^2, 2.756e-6 N/cm^2, 43.373], pmax = 73 kPa = 7.3 N/cm^2
# Heidlauf16: [3.176e-10 N/cm^2, 1.813 N/cm^2, 1.075e-2 N/cm^2, 9.1733], pmax = 7.3 N/cm^2

c1 = 3.176e-10              # [N/cm^2]
c2 = 1.813                  # [N/cm^2]
b  = 1.075e-2               # [N/cm^2] anisotropy parameter
d  = 9.1733                 # [-] anisotropy parameter
material_parameters = [c1, c2, b, d]   # material parameters
pmax = 7.3                  # [N/cm^2] maximum isometric active stress

constant_body_force = (0,0,-9.81e-4)   # [cm/ms^2], gravity constant for the body force
#constant_body_force = (0,0,0)

# traction force, given by command line arguments --yforce and --zforce
bottom_traction = [0.0,variables.yforce,variables.zforce]    # [N]

# boundary conditions (for quadratic elements)
[mx, my, mz] = variables.meshes["3Dmesh_quadratic"]["nPointsGlobal"]
[nx, ny, nz] = variables.meshes["3Dmesh_quadratic"]["nElements"]

fiber_mesh_names = [mesh_name for mesh_name in variables.meshes.keys() if "MeshFiber" in mesh_name]

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
       
# Neumann BC at bottom nodes, traction downwards
elasticity_neumann_bc = [{"element": 0*nx*ny + j*nx + i, "constantVector": bottom_traction, "face": "2-"} for j in range(ny) for i in range(nx)]

# Neumann boundary conditions, specify upward force for top elements, slightly in y-direction
#neumann_bc = [{"element": (nz-1)*nx*ny + j*nx + i, "constantVector": bottom_traction, "face": "2+"} for j in range(ny) for i in range(nx)]

# function to postprocess the output
def postprocess(result):
  result = result[0]
  #print(result)
  
  current_time = result["currentTime"]
  timestep_no = result["timeStepNo"]
  print("timestep_no",timestep_no)
  if timestep_no == 0:
    return
  
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
    f.write("nonlinear;{};{}\n".format(variables.zforce,uz[0]))

config = {
  "scenarioName": "3d_muscle",
  "logFormat":    "csv",     # "csv" or "json", format of the lines in the log file, csv gives smaller files
  "solverStructureDiagramFile":     "out/" + variables.scenario_name + "/solver_structure.txt",     # output file of a diagram that shows data connection between solvers
  "mappingsBetweenMeshesLogFile":   "out/" + variables.scenario_name + "/mappings_between_meshes.txt",   # log file for mappings between meshes
  
  "Meshes": variables.meshes,
  "HyperelasticitySolver": {
    "durationLogKey": "nonlinear",
    
    "materialParameters": material_parameters,
    "displacementsScalingFactor": 1.0,            # scaling factor for output of displacements
    "constantBodyForce": constant_body_force,     # body force in whole body region
    "residualNormLogFilename": "log_residual_norm.txt",
    "useAnalyticJacobian": True,
    "useNumericJacobian": False,   # Only works in parallel execution. If both numeric and analytic are enable, it uses the analytic for the preconditioner and the numeric as normal jacobian
      
    "dumpDenseMatlabVariables": False,   # extra output of matlab vectors, x,r, jacobian matrix
    # if useAnalyticJacobian,useNumericJacobian and dumpDenseMatlabVariables are all three true, the analytic and numeric jacobian matrices will get compared to see if there are programming errors for the analytic jacobian
    
    # mesh
    "meshName": "3Dmesh_quadratic",      # mesh with quadratic Lagrange ansatz functions
    "inputMeshIsGlobal": True,           # boundary conditions are specified in global numberings
    
    "fiberMeshNames": fiber_mesh_names,   # fiber meshes that will be used to determine the fiber direction
    
    # nonlinear solver
    "relativeTolerance": 1e-5,          # 1e-10 relative tolerance of the linear solver
    "absoluteTolerance": 1e-5,          # 1e-10 absolute tolerance of the residual of the linear solver    
    "solverType": "preonly",            # type of the linear solver: cg groppcg pipecg pipecgrr cgne nash stcg gltr richardson chebyshev gmres tcqmr fcg pipefcg bcgs ibcgs fbcgs fbcgsr bcgsl cgs tfqmr cr pipecr lsqr preonly qcg bicg fgmres pipefgmres minres symmlq lgmres lcd gcr pipegcr pgmres dgmres tsirm cgls
    "preconditionerType": "lu",         # type of the preconditioner
    "maxIterations": 1e4,               # maximum number of iterations in the linear solver
    "dumpFilename": "",#"out/m",            # filename for output of solver matrix
    "dumpFormat": "matlab",             # default, ascii, matlab
    "snesMaxFunctionEvaluations": 1e8,  # maximum number of function iterations
    "snesMaxIterations": 15,             # maximum number of iterations in the nonlinear solver
    "snesRebuildJacobianFrequency": variables.njacobi,  # frequency with which the jacobian is newly computed
    "snesRelativeTolerance": 1e-5,      # relative tolerance of the nonlinear solver
    "snesLineSearchType": "l2",         # type of linesearch, possible values: "bt" "nleqerr" "basic" "l2" "cp" "ncglinear"
    "snesAbsoluteTolerance": 1e-5,      # absolute tolerance of the nonlinear solver
    "loadFactorGiveUpThreshold": 0.1,   # if the adaptive time stepping produces a load factor smaller than this value, the solution will be accepted for the current timestep, even if it did not converge fully to the tolerance
    
    #"loadFactors":  [0.1, 0.2, 0.35, 0.5, 1.0],   # load factors for every timestep
    "loadFactors": [],                  # no load factors, solve problem directly
    "scaleInitialGuess":          False,# when load stepping is used, scale initial guess between load steps a and b by sqrt(a*b)/a. This potentially reduces the number of iterations per load step (but not always).
    "nNonlinearSolveCalls": 1,          # how often the nonlinear solve should be repeated
    
    # boundary conditions
    "dirichletBoundaryConditions": elasticity_dirichlet_bc,
    "neumannBoundaryConditions": elasticity_neumann_bc,
    "divideNeumannBoundaryConditionValuesByTotalArea": True,            # if the given Neumann boundary condition values under "neumannBoundaryConditions" are total forces instead of surface loads and therefore should be scaled by the surface area of all elements where Neumann BC are applied
    #"updateDirichletBoundaryConditionsFunction": update_dirichlet_boundary_conditions,
    "updateDirichletBoundaryConditionsFunction": None,
    "updateDirichletBoundaryConditionsFunctionCallInterval": 1,
    
    "dirichletOutputFilename":     "out/" + variables.scenario_name + "/dirichlet_boundary_conditions",    # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
      
    "OutputWriter" : [   # output files for displacements function space (quadratic elements)
      {"format": "Paraview", "outputInterval": 1, "filename": "out/" + variables.scenario_name + "/u", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
      {"format": "PythonCallback", "callback": postprocess, "outputInterval": 1, "binary":True, "onlyNodalValues":True, "fileNumbering": "incremental"},
    ],
    "pressure": {   # output files for pressure function space (linear elements)
      "OutputWriter" : [
        {"format": "Paraview", "outputInterval": 1, "filename": "out/" + variables.scenario_name + "/p", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
        {"format": "PythonFile", "filename": "out/" + variables.scenario_name + "/p", "outputInterval": 1, "binary":False, "onlyNodalValues":True, "fileNumbering": "incremental"},
      ]
    },
    # output writer for debugging, outputs files after each load increment, the geometry is not changed but u and v are written
    "LoadIncrements": {   
      "OutputWriter" : [
        {"format": "Paraview", "outputInterval": 1, "filename": "out/" + variables.scenario_name + "/load_increments", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
      ]
    },
  },
}
