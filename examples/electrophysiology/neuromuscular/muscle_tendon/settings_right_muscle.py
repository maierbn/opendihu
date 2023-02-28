import sys, os
import timeit
import argparse
import importlib
import distutils.util

# parse rank arguments
rank_no = int(sys.argv[-2])
n_ranks = int(sys.argv[-1])

# add variables subfolder to python path where the variables script is located
script_path = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, script_path)
sys.path.insert(0, os.path.join(script_path,'variables'))

import variables              # file variables.py, defined default values for all parameters, you can set the parameters there


variables.n_subdomains = variables.n_subdomains_x*variables.n_subdomains_y*variables.n_subdomains_z



# automatically initialize partitioning if it has not been set
if n_ranks != variables.n_subdomains:

  # create all possible partitionings to the given number of ranks
  optimal_value = n_ranks**(1/3)
  possible_partitionings = []
  for i in range(1,n_ranks+1):
    for j in range(1,n_ranks+1):
      if i*j <= n_ranks and n_ranks % (i*j) == 0:
        k = int(n_ranks / (i*j))
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

# output information of run
if rank_no == 0:
  print("scenario_name: {},  n_subdomains: {} {} {},  n_ranks: {},  end_time: {}".format(variables.scenario_name, variables.n_subdomains_x, variables.n_subdomains_y, variables.n_subdomains_z, n_ranks, variables.end_time))
  print("cellml_file:             {}".format(variables.cellml_file))
  print("firing_times_file:       {}".format(variables.firing_times_file))
  print("********************************************************************************")

  # start timer to measure duration of parsing of this script
  t_start_script = timeit.default_timer()

# initialize all helper variables
from helper import *

variables.scenario_name = "muscle_right"

variables.n_subdomains_xy = variables.n_subdomains_x * variables.n_subdomains_y
variables.n_fibers_total = variables.n_fibers_x * variables.n_fibers_y



#### set Dirichlet BC and Neumann BC for the free side of the muscle

[nx, ny, nz] = [elem + 1 for elem in variables.n_elements_muscle1]
[mx, my, mz] = [elem // 2 for elem in variables.n_elements_muscle1] # quadratic elements consist of 2 linear elements along each axis

k = nz -1  #fixed side of the muscle

for j in range(ny):
    for i in range(nx):
      variables.elasticity_dirichlet_bc[k*nx*ny + j*nx + i] = [0.0, 0.0, 0.0, None,None,None] # displacement ux uy uz, velocity vx vy vz


# meshes

# add neuron meshes
muscle_meshes = {

  "muscle2Mesh": {
    "nElements" :         variables.n_elements_muscle1,
    "physicalExtent":     variables.muscle1_extent,
    "physicalOffset":     [0, 0, variables.muscle1_extent[2] + variables.tendon_extent[2]],
    "logKey":             "muscle2",
    "inputMeshIsGlobal":  True,
    "nRanks":             n_ranks
  },

  # needed for mechanics solver
  "muscle2Mesh_quadratic": {
    "nElements" :         [elems // 2 for elems in variables.n_elements_muscle1],
    "physicalExtent":     variables.muscle1_extent,
    "physicalOffset":     [0, 0, variables.muscle1_extent[2] + variables.tendon_extent[2]],
    "logKey":             "muscle2_quadratic",
    "inputMeshIsGlobal":  True,
    "nRanks":             n_ranks,
  }
}
variables.meshes.update(muscle_meshes)
variables.meshes.update(fiber_meshes)


def dbg(x, name=None):
  if name:
    print(name, end=': ')
  print(x)
  return x

# define the config dict
config = {
  "scenarioName":                   variables.scenario_name,      # scenario name to identify the simulation runs in the log file
  "logFormat":                      "csv",                        # "csv" or "json", format of the lines in the log file, csv gives smaller files
  "solverStructureDiagramFile":     "solver_structure.txt",       # output file of a diagram that shows data connection between solvers
  "mappingsBetweenMeshesLogFile":   "mappings_between_meshes_log.txt",    # log file for mappings 
  "Meshes":                         variables.meshes,

  "PreciceAdapter": {        # precice adapter for bottom tendon
      "timeStepOutputInterval":   100,                        # interval in which to display current timestep and time in console
      "timestepWidth":            variables.dt_elasticity,                          # coupling time step width, must match the value in the precice config
      "couplingEnabled":          True,                       # if the precice coupling is enabled, if not, it simply calls the nested solver, for debugging
      "preciceConfigFilename":    "../precice_config_multi_coupling.xml",    # the preCICE configuration file
      "preciceParticipantName":   "MuscleSolverRight",             # name of the own precice participant, has to match the name given in the precice xml config file
      "preciceMeshes": [                                      # the precice meshes get created as the top or bottom surface of the main geometry mesh of the nested solver
        {
          "preciceMeshName":      "MuscleMeshRight",            # precice name of the 2D coupling mesh
          "face":                 "2-",                       # face of the 3D mesh where the 2D mesh is located, "2-" = bottom, "2+" = top
        }
      ],
      "preciceData": [  
        {
          "mode":                 "write-displacements-velocities",   # mode is one of "read-displacements-velocities", "read-traction", "write-displacements-velocities", "write-traction"
          "preciceMeshName":      "MuscleMeshRight",                    # name of the precice coupling surface mesh, as given in the precice xml settings file
          "displacementsName":    "Displacement",                     # name of the displacements "data", i.e. field variable, as given in the precice xml settings file
          "velocitiesName":       "Velocity",                     # name of the velocities "data", i.e. field variable, as given in the precice xml settings file

        },
        {
          "mode":                 "read-traction",                    # mode is one of "read-displacements-velocities", "read-traction", "write-displacements-velocities", "write-traction"
          "preciceMeshName":      "MuscleMeshRight",                    # name of the precice coupling surface mesh, as given in the precice xml settings 
          "tractionName":         "Traction",                         # name of the traction "data", i.e. field variable, as given in the precice xml settings file
        }
      ],
    
    "DynamicHyperelasticitySolver": {
      "timeStepWidth":              variables.dt_elasticity,#variables.dt_elasticity,      # time step width 
      "endTime":                    variables.end_time,           # end time of the simulation time span    
      "durationLogKey":             "duration_mechanics",         # key to find duration of this solver in the log file
      "timeStepOutputInterval":     1,                            # how often the current time step should be printed to console
      
      "materialParameters":         variables.muscle_material_parameters,  # material parameters of the Mooney-Rivlin material
      "density":                    variables.rho,  
      "displacementsScalingFactor": 1.0,                          # scaling factor for displacements, only set to sth. other than 1 only to increase visual appearance for very small displacements
      "residualNormLogFilename":    "log_residual_norm.txt",      # log file where residual norm values of the nonlinear solver will be written
      "useAnalyticJacobian":        True,                         # whether to use the analytically computed jacobian matrix in the nonlinear solver (fast)
      "useNumericJacobian":         False,                        # whether to use the numerically computed jacobian matrix in the nonlinear solver (slow), only works with non-nested matrices, if both numeric and analytic are enable, it uses the analytic for the preconditioner and the numeric as normal jacobian
        
      "dumpDenseMatlabVariables":   False,                        # whether to have extra output of matlab vectors, x,r, jacobian matrix (very slow)
      # if useAnalyticJacobian,useNumericJacobian and dumpDenseMatlabVariables all all three true, the analytic and numeric jacobian matrices will get compared to see if there are programming errors for the analytic jacobian
      
      # mesh
      "meshName":                   "muscle2Mesh_quadratic",           # mesh with quadratic Lagrange ansatz functions
      "inputMeshIsGlobal":          True,                         # boundary conditions are specified in global numberings, whereas the mesh is given in local numberings
      
      "fiberMeshNames":             [],                           # fiber meshes that will be used to determine the fiber direction
      "fiberDirection":             [0,0,1],                      # if fiberMeshNames is empty, directly set the constant fiber direction, in element coordinate system
      
      # nonlinear solver
      "relativeTolerance":          1e-10,                         # 1e-10 relative tolerance of the linear solver
      "absoluteTolerance":          1e-10,                        # 1e-10 absolute tolerance of the residual of the linear solver       
      "solverType":                 "preonly",                    # type of the linear solver: cg groppcg pipecg pipecgrr cgne nash stcg gltr richardson chebyshev gmres tcqmr fcg pipefcg bcgs ibcgs fbcgs fbcgsr bcgsl cgs tfqmr cr pipecr lsqr preonly qcg bicg fgmres pipefgmres minres symmlq lgmres lcd gcr pipegcr pgmres dgmres tsirm cgls
      "preconditionerType":         "lu",                         # type of the preconditioner
      "maxIterations":              1e4,                          # maximum number of iterations in the linear solver
      "snesMaxFunctionEvaluations": 1e8,                          # maximum number of function iterations
      "snesMaxIterations":          240,                           # maximum number of iterations in the nonlinear solver
      "snesRelativeTolerance":      1e-4,                         # relative tolerance of the nonlinear solver
      "snesLineSearchType":         "l2",                         # type of linesearch, possible values: "bt" "nleqerr" "basic" "l2" "cp" "ncglinear"
      "snesAbsoluteTolerance":      1e-5,                         # absolute tolerance of the nonlinear solver
      "snesRebuildJacobianFrequency": 5,          
      
      #"dumpFilename": "out/r{}/m".format(sys.argv[-1]),          # dump system matrix and right hand side after every solve
      "dumpFilename":               "",                           # dump disabled
      "dumpFormat":                 "matlab",                     # default, ascii, matlab
      
      #"loadFactors":                [0.1, 0.2, 0.35, 0.5, 1.0],   # load factors for every timestep
      #"loadFactors":                [0.5, 1.0],                   # load factors for every timestep
      "loadFactors":                [],                           # no load factors, solve problem directly
      "loadFactorGiveUpThreshold":  1e-3, 
      "nNonlinearSolveCalls":       1,                            # how often the nonlinear solve should be called
      
      # boundary and initial conditions
      "dirichletBoundaryConditions": variables.elasticity_dirichlet_bc,   # the initial Dirichlet boundary conditions that define values for displacements u and velocity v
      "neumannBoundaryConditions":   [],     # Neumann boundary conditions that define traction forces on surfaces of elements
      "divideNeumannBoundaryConditionValuesByTotalArea": False,    # if the initial values for the dynamic nonlinear problem should be computed by extrapolating the previous displacements and velocities
      "updateDirichletBoundaryConditionsFunction": None, #update_dirichlet_bc,   # function that updates the dirichlet BCs while the simulation is running
      "updateDirichletBoundaryConditionsFunctionCallInterval": 1,         # stide every which step the update function should be called, 1 means every time step
      "updateNeumannBoundaryConditionsFunction": None,       # a callback function to periodically update the Neumann boundary conditions
      "updateNeumannBoundaryConditionsFunctionCallInterval": 1,           # every which step the update function should be called, 1 means every time step 
      
      "constantBodyForce":           None,       # a constant force that acts on the whole body, e.g. for gravity
      "initialValuesDisplacements":  [[0.0,0.0,0.0] for _ in range(nx * ny * nz)],     # the initial values for the displacements, vector of values for every node [[node1-x,y,z], [node2-x,y,z], ...]
      "initialValuesVelocities":     [[0.0,0.0,0.0] for _ in range(nx * ny * nz)],     # the initial values for the velocities, vector of values for every node [[node1-x,y,z], [node2-x,y,z], ...]
      "extrapolateInitialGuess":     True, 

      "dirichletOutputFilename":     "out/"+variables.scenario_name+"/dirichlet_boundary_conditions_tendon",    # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
      # "totalForceLogFilename":       "out/tendon_force.csv",              # filename of a log file that will contain the total (bearing) forces and moments at the top and bottom of the volume
      # "totalForceLogOutputInterval": 10,                                  # output interval when to write the totalForceLog file
      # "totalForceBottomElementNosGlobal":  [j*nx + i for j in range(ny) for i in range(nx)],                  # global element nos of the bottom elements used to compute the total forces in the log file totalForceLogFilename
      # "totalForceTopElementNosGlobal":     [(nz-1)*ny*nx + j*nx + i for j in range(ny) for i in range(nx)],   # global element nos of the top elements used to compute the total forces in the log file totalForceTopElementsGlobal

      "OutputWriter" : [
          {"format": "Paraview", "outputInterval": 1, "filename": "out/" + variables.scenario_name + "/mechanics_3D", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
          {"format": "PythonCallback", "outputInterval": 1, "callback": variables.muscle_right_write_to_file, "onlyNodalValues":True, "filename": "", "fileNumbering":'incremental'},

        ],
      # define which file formats should be written
      # 1. main output writer that writes output files using the quadratic elements function space. Writes displacements, velocities and PK2 stresses.
      # 2. additional output writer that writes also the hydrostatic pressure
      "pressure": {   # output files for pressure function space (linear elements), contains pressure values, as well as displacements and velocities
        "OutputWriter" : [
          # {"format": "Paraview", "outputInterval": 1, "filename": "out/"+variables.scenario_name+"/p", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
        ]
      },
      # 3. additional output writer that writes virtual work terms
      "dynamic": {    # output of the dynamic solver, has additional virtual work values 
        "OutputWriter" : [   # output files for displacements function space (quadratic elements)
          # {"format": "Paraview", "outputInterval": 1, "filename": "out/"+variables.scenario_name+"/virtual_work", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
          # {"format": "Paraview", "outputInterval": 1, "filename": "out/"+variables.scenario_name+"/virtual_work", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
        ],
      },
      # 4. output writer for debugging, outputs files after each load increment, the geometry is not changed but u and v are written
      "LoadIncrements": {   
        "OutputWriter" : [
          #{"format": "Paraview", "outputInterval": 1, "filename": "out/load_increments", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
        ]
      },
    }
  }
}

# stop   timer and calculate how long parsing lasted
if rank_no == 0:
  t_stop_script = timeit.default_timer()
  print("Python config parsed in {:.1f}s.".format(t_stop_script - t_start_script))
