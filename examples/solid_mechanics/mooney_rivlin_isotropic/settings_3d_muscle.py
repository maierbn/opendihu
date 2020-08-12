# biceps
#

import numpy as np
import pickle
import sys
sys.path.insert(0, '..')
import variables              # file variables.py, defined default values for all parameters, you can set the parameters there
from create_partitioned_meshes_for_settings import *   # file create_partitioned_meshes_for_settings with helper functions about own subdomain


# input mesh file
fiber_file = "../../../electrophysiology/input/left_biceps_brachii_13x13fibers.bin"

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
sampling_stride_z = 500

sampling_stride_x = 1
sampling_stride_y = 1
sampling_stride_z = 100

# create the partitioning using the script in create_partitioned_meshes_for_settings.py
result = create_partitioned_meshes_for_settings(
    n_subdomains_x, n_subdomains_y, n_subdomains_z, 
    fiber_file, load_fiber_data,
    sampling_stride_x, sampling_stride_y, sampling_stride_z, generate_linear_3d_mesh=True, generate_quadratic_3d_mesh=True)

#parse result
[variables.meshes, variables.own_subdomain_coordinate_x, variables.own_subdomain_coordinate_y, variables.own_subdomain_coordinate_z, variables.n_fibers_x, variables.n_fibers_y, variables.n_points_whole_fiber] = result

node_positions = variables.meshes["3Dmesh_quadratic"]["nodePositions"]
#node_positions = variables.meshes["3Dmesh"]["nodePositions"]

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
material_parameters = [c1, c2]   # material parameters
pmax = 7.3                  # [N/cm^2] maximum isometric active stress

constant_body_force = (0,0,-9.81e-4)   # [cm/ms^2], gravity constant for the body force
#constant_body_force = (0,0,0)
bottom_traction = [0.0,-1,-5]        # [N]
#bottom_traction = [0.0,0.0,0.0]        # [N]


# boundary conditions (for quadratic elements)
[mx, my, mz] = variables.meshes["3Dmesh_quadratic"]["nPointsGlobal"]
[nx, ny, nz] = variables.meshes["3Dmesh_quadratic"]["nElements"]

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


config = {
  "scenarioName": "3d_muscle",
  "logFormat":    "csv",     # "csv" or "json", format of the lines in the log file, csv gives smaller files
  "solverStructureDiagramFile":     "solver_structure.txt",     # output file of a diagram that shows data connection between solvers
  "mappingsBetweenMeshesLogFile":   "mappings_between_meshes.txt",   # log file for mappings between meshes
  
  "Meshes": variables.meshes,
  "HyperelasticitySolver": {
    "durationLogKey": "nonlinear",
    
    #"materialParameters": [1.5,2.0],
    "materialParameters": material_parameters,
    "displacementsScalingFactor": 1.0,            # scaling factor for output of displacements
    "constantBodyForce": constant_body_force,     # body force in whole body region
    "residualNormLogFilename": "log_residual_norm.txt",
    "useAnalyticJacobian": True,
    "useNumericJacobian": False,   # Only works in parallel execution. If both numeric and analytic are enable, it uses the analytic for the preconditioner and the numeric as normal jacobian
      
    "dumpDenseMatlabVariables": False,   # extra output of matlab vectors, x,r, jacobian matrix
    # if useAnalyticJacobian,useNumericJacobian and dumpDenseMatlabVariables are all three true, the analytic and numeric jacobian matrices will get compared to see if there are programming errors for the analytic jacobian
    
    # mesh
    "meshName": "3Dmesh_quadratic",     # mesh with quadratic Lagrange ansatz functions
    "inputMeshIsGlobal": True,          # boundary conditions are specified in global numberings
    
    # nonlinear solver
    "relativeTolerance": 1e-5,         # 1e-10 relative tolerance of the linear solver
    "absoluteTolerance": 1e-5,         # 1e-10 absolute tolerance of the residual of the linear solver    
    "solverType": "preonly",            # type of the linear solver: cg groppcg pipecg pipecgrr cgne nash stcg gltr richardson chebyshev gmres tcqmr fcg pipefcg bcgs ibcgs fbcgs fbcgsr bcgsl cgs tfqmr cr pipecr lsqr preonly qcg bicg fgmres pipefgmres minres symmlq lgmres lcd gcr pipegcr pgmres dgmres tsirm cgls
    "preconditionerType": "lu",         # type of the preconditioner
    "maxIterations": 1e4,               # maximum number of iterations in the linear solver
    "dumpFilename": "",#"out/m",            # filename for output of solver matrix
    "dumpFormat": "matlab",             # default, ascii, matlab
    "snesMaxFunctionEvaluations": 1e8,  # maximum number of function iterations
    "snesMaxIterations": 15,             # maximum number of iterations in the nonlinear solver
    "snesRebuildJacobianFrequency": 5,  # frequency with which the jacobian is newly computed
    "snesRelativeTolerance": 1e-5,      # relative tolerance of the nonlinear solver
    "snesLineSearchType": "l2",         # type of linesearch, possible values: "bt" "nleqerr" "basic" "l2" "cp" "ncglinear"
    "snesAbsoluteTolerance": 1e-5,      # absolute tolerance of the nonlinear solver
    "loadFactorGiveUpThreshold": 0.1,   # if the adaptive time stepping produces a load factor smaller than this value, the solution will be accepted for the current timestep, even if it did not converge fully to the tolerance
    
    #"loadFactors":  [0.1, 0.2, 0.35, 0.5, 1.0],   # load factors for every timestep
    "loadFactors": [],                 # no load factors, solve problem directly
    "nNonlinearSolveCalls": 1,         # how often the nonlinear solve should be repeated
    
    # boundary conditions
    "dirichletBoundaryConditions": elasticity_dirichlet_bc,
    "neumannBoundaryConditions": elasticity_neumann_bc,
    "divideNeumannBoundaryConditionValuesByTotalArea": True,            # if the given Neumann boundary condition values under "neumannBoundaryConditions" are total forces instead of surface loads and therefore should be scaled by the surface area of all elements where Neumann BC are applied
    #"updateDirichletBoundaryConditionsFunction": update_dirichlet_boundary_conditions,
    "updateDirichletBoundaryConditionsFunction": None,
    "updateDirichletBoundaryConditionsFunctionCallInterval": 1,
    
    "dirichletOutputFilename":     "out/dirichlet_boundary_conditions",                                # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
    
    "OutputWriter" : [   # output files for displacements function space (quadratic elements)
      {"format": "Paraview", "outputInterval": 1, "filename": "out/u", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
      {"format": "PythonFile", "filename": "out/u", "outputInterval": 1, "binary":False, "onlyNodalValues":True, "fileNumbering": "incremental"},
    ],
    "pressure": {   # output files for pressure function space (linear elements)
      "OutputWriter" : [
        {"format": "Paraview", "outputInterval": 1, "filename": "out/p", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
        {"format": "PythonFile", "filename": "out/p", "outputInterval": 1, "binary":False, "onlyNodalValues":True, "fileNumbering": "incremental"},
      ]
    },
    # output writer for debugging, outputs files after each load increment, the geometry is not changed but u and v are written
    "LoadIncrements": {   
      "OutputWriter" : [
        {"format": "Paraview", "outputInterval": 1, "filename": "out/load_increments", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
      ]
    },
  },
}
