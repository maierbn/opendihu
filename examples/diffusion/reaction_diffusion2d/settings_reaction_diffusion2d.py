# Diffusion 2D
import numpy as np
import scipy.stats

n = 4   # number of elements

dt = 1e-3   # timestep width

# initial values
iv = {}

for y in range(int(0.2*n), int(0.3*n)):
  for x in range(int(0.5*n), int(0.8*n)):
    i = y*(2*n+1) + x
    iv[i] = 1.0

def set_reaction_term_values(n_dofs_global, n_nodes_global_per_coordinate_direction, time_step_no, current_time, values, global_natural_dofs, additional_argument):
  # n_dofs_global:       (int) global number of dofs in the mesh where to set the values
  # n_nodes_global_per_coordinate_direction (list of ints)   [mx, my, mz] number of global nodes in each coordinate direction. 
  #                       For composite meshes, the values are only for the first submesh, for other meshes sum(...) equals n_dofs_global
  # time_step_no:        (int)   current time step number
  # current_time:        (float) the current simulation time
  # values:              (list of floats) all current local values of the field variable, if there are multiple components, they are stored in struct-of-array memory layout 
  #                       i.e. [point0_component0, point0_component1, ... pointN_component0, point1_component0, point1_component1, ...]
  # global_natural_dofs  (list of ints) for every local dof no. the dof no. in global natural ordering
  # additional_argument: The value of the option "additionalArgument", can be any Python object.
  
  for local_dof_no in range(len(values)):
    global_dof_no = global_natural_dofs[local_dof_no]
    
    mx = n_nodes_global_per_coordinate_direction[0]
    
    i = global_dof_no % mx        # index in x direction
    j = int(global_dof_no / mx)   # index in y direction
    
    # set gaussian
    center = np.array((0.7*mx,0.7*mx))
    x = np.linalg.norm(np.array((i,j)) - center)
    values[local_dof_no] += scipy.stats.norm.pdf(x, scale=0.1*mx)*0.001

config = {
  "solverStructureDiagramFile":     "solver_structure.txt",     # output file of a diagram that shows data connection between solvers
  "logFormat":                      "csv",                      # "csv" or "json", format of the lines in the log file, csv gives smaller files
  "scenarioName":                   "diffusion",                # scenario name to find the run in the log file
  "mappingsBetweenMeshesLogFile":   "out/mappings_between_meshes.txt",  # log file for mappings between meshes
  "Meshes": {
    "mesh": {
      "nElements":         [n,n],         # number of elements in x and y direction
      "physicalExtent":    [4.0,4.0],     # the size of the domain in physical space
      "physicalOffset":    [0.0,0.0],     # position of the lower left point of the mesh
      "inputMeshIsGlobal": True,          # if the values in "nElements" and "physicalExtent" are global numbers
    }
  },
  "Solvers": {
    "solver": {
      "relativeTolerance":  1e-10,        # relative tolerance of the residual normal, respective to the initial residual norm, linear solver
      "absoluteTolerance":  1e-10,        # 1e-10 absolute tolerance of the residual          
      "maxIterations":      1e5,          # maximum number of iterations
      "solverType":         "gmres",      # type of the linear solver
      "preconditionerType": "none",       # type of the preconditioner
      "dumpFormat":         "default",    # one of "ascii", "matlab", "default", format if matrix and rhs should be written to a file after each solve
      "dumpFilename":       "",           # filename if matrix and rhs should be written to a file after each solve, if empty, no file will be written
    },
  },
  "GodunovSplitting": {
    "timeStepWidth":          dt,             # time step width of the coupling, i.e. when the diffusion and reaction terms will be evaluated
    "logTimeStepWidthAsKey":  "dt_coupling",       # the key under which the timeStepWidth value will be written to the log file
    "durationLogKey":         "duration_total",    # a key to find the total computational duration in the log file
    "timeStepOutputInterval": 10,             # how often to print information about iteration to the console
    "endTime":                10,             # end time of the simulation
    
    "connectedSlotsTerm1To2": {0:0},          # diffusion term solution <=> reaction term solution
    "connectedSlotsTerm2To1": {0:0},          # diffusion term solution <=> reaction term solution

    "Term1": {      # diffusion term
      "ExplicitEuler": {
        "initialValues":      iv,             # initial values to be set
        "inputMeshIsGlobal":  True,           # if the initial values are global numbers
        "timeStepWidth":      dt,             # potentially subcycling for diffusion term, here we set the diffusion timestep width and the outer, coupling timestep width to be the same
        "endTime":            1.0,            # the end time value here will be ignored and instead set by the Godunov scheme to match the outer time step
        "checkForNanInf":     False,          # if the solution should be checked for nan and inf values, after each iteration.
        "timeStepOutputInterval": 1,          # how often to print information about iteration to the console
        "dirichletBoundaryConditions": {},    # dirichlet boundary conditions
        "dirichletOutputFilename": None,      # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
        "nAdditionalFieldVariables":  0,      # how many additional field variables there should be that can be connected to output slots and that will be written to the output file
        "additionalSlotNames":  [],           # slot names for the additional field variables
        
        "FiniteElementMethod" : {
          "meshName":          "mesh",        # reference to the mesh, defined under "Meshes"
          "solverName":        "solver",      # reference to the linear system solver
          "prefactor":         0.1,           # prefactor c
          "inputMeshIsGlobal": True,          # if the boundary conditions are given as global numbers
          "slotName":          "",            # slot name, needed if slots are connected by the global option "connectedSlots" and thereby identified via slot name
        },
        "OutputWriter" : [
          #{"format": "Paraview", "interval": 1, "filename": "out", "binaryOutput": "false", "fixedFormat": False, "onlyNodalValues":True, "fileNumbering": "incremental"},
          {"format": "PythonFile", "outputInterval": 100, "filename": "out/diffusion2d", "onlyNodalValues":True, "binary": True, "fileNumbering": "incremental"}
        ]
      },
    },
    "Term2": {     # reaction term
      "PrescribedValues": {
        "meshName":            "mesh",        # the mesh to use for the field variables
        "numberTimeSteps":     1,             # number of timesteps to call the callback functions subsequently, this is usually 1 for prescribed values, because it is enough to set the reaction term only once per time step
        "timeStepOutputInterval": 10,         # if the time step should be written to console, a value > 1 produces no output
        "slotNames":           [],            # slot names for the connector slots of the set field variables, optional
        
        # a list of field variables that will get values assigned in every timestep, by the provided callback function
        "fieldVariables1": [
          {"name": "reaction_term", "callback": set_reaction_term_values},
        ],
        "fieldVariables2":     [],
        "additionalArgument":  None,          # a custom argument to the fieldVariables callback functions, this will be passed on as the last argument
        
        "OutputWriter" : [
          #{"format": "Paraview", "interval": 1, "filename": "out", "binaryOutput": "false", "fixedFormat": False, "onlyNodalValues":True, "fileNumbering": "incremental"},
          {"format": "PythonFile", "outputInterval": 100, "filename": "out/reaction", "binary": True, "onlyNodalValues":True, "fileNumbering": "incremental"}
        ]
      },
    },
  }
}
