# Transversely-isotropic Mooney Rivlin on a tendon geometry
# Note, this is not possible to be run in parallel because the fibers cannot be initialized without MultipleInstances class.
import sys, os
import numpy as np
import sys

# parse rank arguments
rank_no = (int)(sys.argv[-2])
n_ranks = (int)(sys.argv[-1])

#add variables subfolder to python path where the variables script is located
script_path = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, script_path)
sys.path.insert(0, os.path.join(script_path,'variables'))

import variables              # file variables.py, defines default values for all parameters, you can set the parameters there
from create_partitioned_meshes_for_settings import *   # file create_partitioned_meshes_for_settings with helper functions about own subdomain
from helper import *

# modify variables according to specific scenario

variables.scenario_name = "tendon"

# material parameters for tendon material
c = 9.98                    # [N/cm^2=kPa]
ca = 14.92                  # [-]
ct = 14.7                   # [-]
cat = 9.64                  # [-]
ctt = 11.24                 # [-]
mu = 3.76                   # [N/cm^2=kPa]
k1 = 42.217e3               # [N/cm^2=kPa]
k2 = 411.360e3              # [N/cm^2=kPa]
variables.material_parameters = [c, ca, ct, cat, ctt, mu, k1, k2]

 

# compute partitioning
if rank_no == 0:
  if n_ranks != variables.n_subdomains_x*variables.n_subdomains_y*variables.n_subdomains_z:
    print("\n\nError! Number of ranks {} does not match given partitioning {} x {} x {} = {}.\n\n".format(n_ranks, variables.n_subdomains_x, variables.n_subdomains_y, variables.n_subdomains_z, variables.n_subdomains_x*variables.n_subdomains_y*variables.n_subdomains_z))
    sys.exit(-1)
    

# add meshes
meshes_tendon = {
  # no `nodePositions` fields as the nodes are created internally
  "tendon_Mesh": {
    "nElements" :         variables.n_elements_tendon,
    "physicalExtent":     variables.tendon_extent,
    "physicalOffset":     variables.tendon_offset,
    "logKey":             "tendon",
    "inputMeshIsGlobal":  True,
    "nRanks":             n_ranks
  },
  # needed for mechanics solver
  "tendon_Mesh_quadratic": {
    "nElements" :         [elems // 2 for elems in variables.n_elements_tendon],
    "physicalExtent":     variables.tendon_extent,
    "physicalOffset":     variables.tendon_offset,
    "logKey":             "tendon_quadratic",
    "inputMeshIsGlobal":  True,
    "nRanks":             n_ranks,
  }
}
variables.meshes.update(meshes_tendon)

# boundary conditions (for quadratic elements)
# --------------------------------------------

[nx, ny, nz] = [elem + 1 for elem in variables.n_elements_tendon]
[mx, my, mz] = [elem // 2 for elem in variables.n_elements_tendon] # quadratic elements consist of 2 linear elements along each axis

variables.elasticity_dirichlet_bc = {}

k = nz-1 #free side of the tendon

for j in range(ny):
    for i in range(nx):
      variables.elasticity_dirichlet_bc[k*nx*ny + j*nx + i] = [None, None, 0.0, None, None, None] # displacement ux uy uz, velocity vx vy vz

for k in range(nz):
    variables.elasticity_dirichlet_bc[k*nx*ny] = [0.0, 0.0, None, None, None, None] # displacement ux uy uz, velocity vx vy vz

variables.elasticity_dirichlet_bc[(nz-1)*nx*ny] = [0.0, 0.0, 0.0, None, None, None] # displacement ux uy uz, velocity vx vy vz


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
    "preciceConfigFilename":    "../precice_config.xml",    # the preCICE configuration file
    "preciceParticipantName":   "TendonSolver",             # name of the own precice participant, has to match the name given in the precice xml config file
    "preciceMeshes": [                                      # the precice meshes get created as the top or bottom surface of the main geometry mesh of the nested solver
      {
        "preciceMeshName":      "TendonMeshLeft",            # precice name of the 2D coupling mesh
        "face":                 "2-",                       # face of the 3D mesh where the 2D mesh is located, "2-" = bottom, "2+" = top
      }
    ],
    "preciceData": [  
      {
        "mode":                 "read-displacements-velocities",   # mode is one of "read-displacements-velocities", "read-traction", "write-displacements-velocities", "write-traction"
        "preciceMeshName":      "TendonMeshLeft",                    # name of the precice coupling surface mesh, as given in the precice xml settings file
        "displacementsName":    "Displacement",                     # name of the displacements "data", i.e. field variable, as given in the precice xml settings file
        "velocitiesName":       "Velocity",                         # name of the velocity "data", i.e. field variable, as given in the precice xml settings file
      },
      {
        "mode":                 "write-traction",                    # mode is one of "read-displacements-velocities", "read-traction", "write-displacements-velocities", "write-traction"
        "preciceMeshName":      "TendonMeshLeft",                    # name of the precice coupling surface mesh, as given in the precice xml settings 
        "tractionName":         "Traction",                         # name of the traction "data", i.e. field variable, as given in the precice xml settings file
      }
    ],
    
    "DynamicHyperelasticitySolver": {
      "timeStepWidth":              variables.dt_elasticity,#variables.dt_elasticity,      # time step width 
      "endTime":                    variables.end_time,           # end time of the simulation time span    
      "durationLogKey":             "duration_mechanics",         # key to find duration of this solver in the log file
      "timeStepOutputInterval":     1,                            # how often the current time step should be printed to console
      
      "materialParameters":         variables.material_parameters,  # material parameters of the Mooney-Rivlin material
      "density":                    variables.rho,  
      "dampingFactor":              variables.damping_factor,  # factor for velocity dependent damping
      "displacementsScalingFactor": 1.0,                          # scaling factor for displacements, only set to sth. other than 1 only to increase visual appearance for very small displacements
      "residualNormLogFilename":    "log_residual_norm.txt",      # log file where residual norm values of the nonlinear solver will be written
      "useAnalyticJacobian":        True,                         # whether to use the analytically computed jacobian matrix in the nonlinear solver (fast)
      "useNumericJacobian":         False,                        # whether to use the numerically computed jacobian matrix in the nonlinear solver (slow), only works with non-nested matrices, if both numeric and analytic are enable, it uses the analytic for the preconditioner and the numeric as normal jacobian
        
      "dumpDenseMatlabVariables":   False,                        # whether to have extra output of matlab vectors, x,r, jacobian matrix (very slow)
      # if useAnalyticJacobian,useNumericJacobian and dumpDenseMatlabVariables all all three true, the analytic and numeric jacobian matrices will get compared to see if there are programming errors for the analytic jacobian
      
      # mesh
      "meshName":                   "tendon_Mesh_quadratic",           # mesh with quadratic Lagrange ansatz functions
      "inputMeshIsGlobal":          True,                         # boundary conditions are specified in global numberings, whereas the mesh is given in local numberings
      
    #   "fiberMeshNames":             [],                           # fiber meshes that will be used to determine the fiber direction
    #   "fiberDirection":             [0,0,1],                      # if fiberMeshNames is empty, directly set the constant fiber direction, in element coordinate system
      
      # nonlinear solver
         # how often the jacobian should be recomputed, -1 indicates NEVER rebuild, 1 means rebuild every time the Jacobian is computed within a single nonlinear solve, 2 means every second time the Jacobian is built etc. -2 means rebuild at next chance but then never again 
      
      #"dumpFilename": "out/r{}/m".format(sys.argv[-1]),          # dump system matrix and right hand side after every solve
      "dumpFilename":               "",                           # dump disabled
      "dumpFormat":                 "matlab",                     # default, ascii, matlab
      
      #"loadFactors":                [0.1, 0.2, 0.35, 0.5, 1.0],   # load factors for every timestep
      #"loadFactors":                [0.5, 1.0],                   # load factors for every timestep
      "loadFactors":                [],                           # no load factors, solve problem directly
      "loadFactorGiveUpThreshold":  1e-5, 
      "scaleInitialGuess":          False,                     # when load stepping is used, scale initial guess between load steps a and b by sqrt(a*b)/a. This potentially reduces the number of iterations per load step (but not always).
      "nNonlinearSolveCalls":       1,                            # how often the nonlinear solve should be called
      
      # boundary and initial conditions
      "dirichletBoundaryConditions": variables.elasticity_dirichlet_bc,   # the initial Dirichlet boundary conditions that define values for displacements u and velocity v
      "neumannBoundaryConditions":   variables.elasticity_neumann_bc,     # Neumann boundary conditions that define traction forces on surfaces of elements
      "divideNeumannBoundaryConditionValuesByTotalArea": False,            # if the given Neumann boundary condition values under "neumannBoundaryConditions" are total forces instead of surface loads and therefore should be scaled by the surface area of all elements where Neumann BC are applied
      "updateDirichletBoundaryConditionsFunction": None, #update_dirichlet_bc,   # function that updates the dirichlet BCs while the simulation is running
      "updateDirichletBoundaryConditionsFunctionCallInterval": 1,         # stide every which step the update function should be called, 1 means every time step
      "updateNeumannBoundaryConditionsFunction": None,       # a callback function to periodically update the Neumann boundary conditions
      "updateNeumannBoundaryConditionsFunctionCallInterval": 1,           # every which step the update function should be called, 1 means every time step      # every which step the update function should be called, 1 means every time step
      
      "initialValuesDisplacements":  [[0.0,0.0,0.0] for _ in range(nx*ny*nz)],     # the initial values for the displacements, vector of values for every node [[node1-x,y,z], [node2-x,y,z], ...]
      "initialValuesVelocities":     [[0.0,0.0,0.0] for _ in range(nx*ny*nz)],     # the initial values for the velocities, vector of values for every node [[node1-x,y,z], [node2-x,y,z], ...]
      "extrapolateInitialGuess":     True,                                # if the initial values for the dynamic nonlinear problem should be computed by extrapolating the previous displacements and velocities
      "constantBodyForce":           variables.constant_body_force,       # a constant force that acts on the whole body, e.g. for gravity
      
      "dirichletOutputFilename":     "out/"+variables.scenario_name+"/dirichlet_boundary_conditions_tendon",    # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
      # "totalForceLogFilename":       "out/tendon_force.csv",              # filename of a log file that will contain the total (bearing) forces and moments at the top and bottom of the volume
      # "totalForceLogOutputInterval": 10,                                  # output interval when to write the totalForceLog file
      # "totalForceBottomElementNosGlobal":  [j*nx + i for j in range(ny) for i in range(nx)],                  # global element nos of the bottom elements used to compute the total forces in the log file totalForceLogFilename
      # "totalForceTopElementNosGlobal":     [(nz-1)*ny*nx + j*nx + i for j in range(ny) for i in range(nx)],   # global element nos of the top elements used to compute the total forces in the log file totalForceTopElementsGlobal

      "OutputWriter" : [
            {"format": "Paraview", "outputInterval": 1, "filename": "out/" + variables.scenario_name + "/mechanics_3D", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
          ],
      # define which file formats should be written
      # 1. main output writer that writes output files using the quadratic elements function space. Writes displacements, velocities and PK2 stresses.
      # 2. additional output writer that writes also the hydrostatic pressure
      "pressure": {   # output files for pressure function space (linear elements), contains pressure values, as well as displacements and velocities
        "OutputWriter" : [
          #{"format": "Paraview", "outputInterval": 1, "filename": "out/"+variables.scenario_name+"/p", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
        ]
      },
      # 3. additional output writer that writes virtual work terms
      "dynamic": {    # output of the dynamic solver, has additional virtual work values 
        "OutputWriter" : [   # output files for displacements function space (quadratic elements)
                {"format": "Paraview", "outputInterval": 1, "filename": "out/"+variables.scenario_name+"/virtual_work", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
          #{"format": "Paraview", "outputInterval": 1, "filename": "out/"+variables.scenario_name+"/virtual_work", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
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