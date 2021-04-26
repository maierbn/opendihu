# This settings file can be used for two different equations:
# - Isotropic hyperelastic material
# - Linear elasticity
#
# arguments: <scenario_name> <force>


import numpy as np
import sys, os

# parameters
force = 1.0                       # [N] load on top
material_parameters = [1.0, 1.0]
physical_extent = [1.0, 1.0, 1.0]
constant_body_force = None
scenario_name = "tensile_test"
dirichlet_bc_mode = "fix_floating"
 
if len(sys.argv) > 3:
  scenario_name = sys.argv[0]
  force = float(sys.argv[1])
  print("scenario_name: {}".format(scenario_name))
  print("force: {}".format(force))
    
  # set material parameters depending on scenario name
  if scenario_name == "compressible_mooney_rivlin":
    material_parameters = [5, 7.5, 10]      # c1, c2, c
    
  elif scenario_name == "compressible_mooney_rivlin_decoupled":
    material_parameters = [5, 7.5, 10.0]      # c1, c2, kappa
    
  elif scenario_name == "incompressible_mooney_rivlin":
    material_parameters = [10, 10]      # c1, c2
    
  elif scenario_name == "nearly_incompressible_mooney_rivlin":
    material_parameters = [10, 10, 1e3]      # c1, c2, kappa

  elif scenario_name == "nearly_incompressible_mooney_rivlin_decoupled":
    material_parameters = [10, 10, 1e3]      # c1, c2, kappa

  elif scenario_name == "linear":
    pass

  elif scenario_name == "nearly_incompressible_mooney_rivlin_febio":
    material_parameters = [10, 10, 1e3]      # c1, c2, kappa

  else:
    print("Error! Please specify the correct scenario, see settings.py for allowed values.\n")
    quit()
 
# number of elements (2x2x2)
nx = 8
ny = 8
nz = 8

# number of nodes
mx = 2*nx + 1
my = 2*ny + 1
mz = 2*nz + 1

# boundary conditions (for quadratic elements)
# --------------------------------------------

# set Dirichlet BC, fix bottom
elasticity_dirichlet_bc = {}
k = 0

# fix z value on the whole x-y-plane
for j in range(my):
  for i in range(mx):
    elasticity_dirichlet_bc[k*mx*my + j*mx + i] = [None,None,0.0,None,None,None]

# fix left edge 
for j in range(my):
  elasticity_dirichlet_bc[k*mx*my + j*mx + 0][0] = 0.0
  
# fix front edge 
for i in range(mx):
  elasticity_dirichlet_bc[k*mx*my + 0*mx + i][1] = 0.0
       
# set Neumann BC, set traction at the top
k = nz-1
traction_vector = [0, 0, force]     # the traction force in specified in the reference configuration

elasticity_neumann_bc = [{"element": k*nx*ny + j*nx + i, "constantVector": traction_vector, "face": "2+"} for j in range(ny) for i in range(nx)]

# callback for result
def handle_result_hyperelasticity(result):
  data = result[0]
  
  if data["timeStepNo"] == 1:
    field_variables = data["data"]
    
    # field_variables[0]: geometry
    # field_variables[1]: u
    # field_variables[2]: v
    # field_variables[3]: t (current traction)
    # field_variables[4]: T (material traction)
    # field_variables[5]: PK2-Stress (Voigt), components: S_11, S_22, S_33, S_12, S_13, S_23
    
    strain = max(field_variables[1]["components"][2]["values"])
    stress = max(field_variables[5]["components"][2]["values"])
    
    print("strain: {}, stress: {}".format(strain, stress))
    
    with open("result.csv","a") as f:
      f.write("{},{},{}\n".format(scenario_name,strain,stress))

# callback for result
def handle_result_febio(result):
  data = result[0]
  
  if data["timeStepNo"] == 1:
    field_variables = data["data"]
    
    strain = max(field_variables[2]["components"][2]["values"])
    stress = max(field_variables[5]["components"][2]["values"])
    
    print("strain: {}, stress: {}".format(strain, stress))
    
    with open("result.csv","a") as f:
      f.write("{},{},{}\n".format(scenario_name,strain,stress))

# callback for result
def handle_result_linear_elasticity(result):
  data = result[0]
  
  if data["timeStepNo"] == -1:
    field_variables = data["data"]
    
    # field_variables[0]: geometry
    # field_variables[1]: solution (displacements)
    # field_variables[2]: rightHandSide
    # field_variables[3]: -rhsNeumannBC
    
    # σ = CC : ε with CC_abcd = K δ_ab δ_cd + μ(δ_ac δ_bd + δ_ad δ_bc - 2/3 δ_ab δ_cd)
    # σ_ab = K*δ_ab*ε_cc + 2*μ*(ε_ab - 1/3*δ_ab*ε_cc)
    # σ_33 = K tr(ε) + μ (ε_33 + ε_33 - 2/3 tr(ε)) =(tensile test in z direction)= (K + 4/3 μ) ε_33
    
    strain = max(field_variables[1]["components"][2]["values"])
    K = 50    # parameters as given in config
    mu = 100
    stress = (K + 4./3*mu) * strain
    
    print("strain: {}, stress: {}".format(strain, stress))
    
    with open("result.csv","a") as f:
      f.write("{},{},{}\n".format(scenario_name,strain,stress))

config = {
  "scenarioName":                 scenario_name,                # scenario name to identify the simulation runs in the log file
  "logFormat":                    "csv",                        # "csv" or "json", format of the lines in the log file, csv gives smaller files
  "solverStructureDiagramFile":   "solver_structure.txt",       # output file of a diagram that shows data connection between solvers
  "mappingsBetweenMeshesLogFile": "mappings_between_meshes_log.txt",    # log file for mappings 
  "Meshes": {
    "3Dmesh_quadratic": { 
      "inputMeshIsGlobal":          True,                       # boundary conditions are specified in global numberings, whereas the mesh is given in local numberings
      "nElements":                  [nx, ny, nz],               # number of quadratic elements in x, y and z direction
      "physicalExtent":             physical_extent,            # physical size of the box
      "physicalOffset":             [0, 0, 0],                  # offset/translation where the whole mesh begins
    },
    "3Dmesh_febio": { 
      "inputMeshIsGlobal":          True,                       # boundary conditions are specified in global numberings, whereas the mesh is given in local numberings
      "nElements":                  [2*nx, 2*ny, 2*nz],               # number of quadratic elements in x, y and z direction
      "physicalExtent":             physical_extent,            # physical size of the box
      "physicalOffset":             [0, 0, 0],                  # offset/translation where the whole mesh begins
    }
  },
  "Solvers": {
    "linearElasticitySolver": {           # solver for linear elasticity
      "relativeTolerance":  1e-10,
      "absoluteTolerance":  1e-10,         # 1e-10 absolute tolerance of the residual    ,
      "maxIterations":      1e4,
      "solverType":         "gmres",
      "preconditionerType": "none",
      "dumpFilename":       "",
      "dumpFormat":         "matlab",
    }, 
  },
  "HyperelasticitySolver": {
    "durationLogKey":             "duration_mechanics",         # key to find duration of this solver in the log file
    
    "materialParameters":         material_parameters,          # material parameters of the Mooney-Rivlin material
    "displacementsScalingFactor": 1.0,                          # scaling factor for displacements, only set to sth. other than 1 only to increase visual appearance for very small displacements
    "residualNormLogFilename":    "log_residual_norm.txt",      # log file where residual norm values of the nonlinear solver will be written
    "useAnalyticJacobian":        True,                         # whether to use the analytically computed jacobian matrix in the nonlinear solver (fast)
    "useNumericJacobian":         False,                        # whether to use the numerically computed jacobian matrix in the nonlinear solver (slow), only works with non-nested matrices, if both numeric and analytic are enable, it uses the analytic for the preconditioner and the numeric as normal jacobian
      
    "dumpDenseMatlabVariables":   False,                        # whether to have extra output of matlab vectors, x,r, jacobian matrix (very slow)
    # if useAnalyticJacobian,useNumericJacobian and dumpDenseMatlabVariables all all three true, the analytic and numeric jacobian matrices will get compared to see if there are programming errors for the analytic jacobian
    
    # mesh
    "meshName":                   "3Dmesh_quadratic",           # mesh with quadratic Lagrange ansatz functions
    "inputMeshIsGlobal":          True,                         # boundary conditions are specified in global numberings, whereas the mesh is given in local numberings
    
    #"fiberMeshNames":             [],                           # fiber meshes that will be used to determine the fiber direction
    #"fiberDirection":             [0,0,1],                      # if fiberMeshNames is empty, directly set the constant fiber direction, in element coordinate system
    
    # nonlinear solver
    "relativeTolerance":          1e-5,                         # 1e-10 relative tolerance of the linear solver
    "absoluteTolerance":          1e-10,                        # 1e-10 absolute tolerance of the residual of the linear solver       
    "solverType":                 "preonly",                    # type of the linear solver: cg groppcg pipecg pipecgrr cgne nash stcg gltr richardson chebyshev gmres tcqmr fcg pipefcg bcgs ibcgs fbcgs fbcgsr bcgsl cgs tfqmr cr pipecr lsqr preonly qcg bicg fgmres pipefgmres minres symmlq lgmres lcd gcr pipegcr pgmres dgmres tsirm cgls
    "preconditionerType":         "lu",                         # type of the preconditioner
    "maxIterations":              1e4,                          # maximum number of iterations in the linear solver
    "snesMaxFunctionEvaluations": 1e8,                          # maximum number of function iterations
    "snesMaxIterations":          100,                           # maximum number of iterations in the nonlinear solver
    "snesRelativeTolerance":      1e-5,                         # relative tolerance of the nonlinear solver
    "snesLineSearchType":         "l2",                         # type of linesearch, possible values: "bt" "nleqerr" "basic" "l2" "cp" "ncglinear"
    "snesAbsoluteTolerance":      1e-5,                         # absolute tolerance of the nonlinear solver
    "snesRebuildJacobianFrequency": 1,                          # how often the jacobian should be recomputed, -1 indicates NEVER rebuild, 1 means rebuild every time the Jacobian is computed within a single nonlinear solve, 2 means every second time the Jacobian is built etc. -2 means rebuild at next chance but then never again 
    
    #"dumpFilename": "out/r{}/m".format(sys.argv[-1]),          # dump system matrix and right hand side after every solve
    "dumpFilename":               "",                           # dump disabled
    "dumpFormat":                 "default",                     # default, ascii, matlab
    
    #"loadFactors":                [0.1, 0.2, 0.35, 0.5, 1.0],   # load factors for every timestep
    #"loadFactors":                [0.5, 1.0],                   # load factors for every timestep
    "loadFactors":                [],                           # no load factors, solve problem directly
    "loadFactorGiveUpThreshold":    0.1,                        # if the adaptive time stepping produces a load factor smaller than this value, the solution will be accepted for the current timestep, even if it did not converge fully to the tolerance
    "nNonlinearSolveCalls":       1,                            # how often the nonlinear solve should be called
    
    # boundary and initial conditions
    "dirichletBoundaryConditions": elasticity_dirichlet_bc,             # the initial Dirichlet boundary conditions that define values for displacements u
    "dirichletOutputFilename":     None,                                # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
    "neumannBoundaryConditions":   elasticity_neumann_bc,               # Neumann boundary conditions that define traction forces on surfaces of elements
    "divideNeumannBoundaryConditionValuesByTotalArea": True,            # if the given Neumann boundary condition values under "neumannBoundaryConditions" are total forces instead of surface loads and therefore should be scaled by the surface area of all elements where Neumann BC are applied
    "updateDirichletBoundaryConditionsFunction": None,                  # function that updates the dirichlet BCs while the simulation is running
    "updateDirichletBoundaryConditionsFunctionCallInterval": 1,         # every which step the update function should be called, 1 means every time step
    
    "initialValuesDisplacements":  [[0.0,0.0,0.0] for _ in range(mx*my*mz)],     # the initial values for the displacements, vector of values for every node [[node1-x,y,z], [node2-x,y,z], ...]
    "initialValuesVelocities":     [[0.0,0.0,0.0] for _ in range(mx*my*mz)],     # the initial values for the velocities, vector of values for every node [[node1-x,y,z], [node2-x,y,z], ...]
    "extrapolateInitialGuess":     True,                                # if the initial values for the dynamic nonlinear problem should be computed by extrapolating the previous displacements and velocities
    "constantBodyForce":           constant_body_force,                 # a constant force that acts on the whole body, e.g. for gravity
    
    "dirichletOutputFilename":      "out/"+scenario_name+"/dirichlet_boundary_conditions",           # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
    
    # define which file formats should be written
    # 1. main output writer that writes output files using the quadratic elements function space. Writes displacements, velocities and PK2 stresses.
    "OutputWriter" : [
      
      # Paraview files
      {"format": "Paraview", "outputInterval": 1, "filename": "out/"+scenario_name+"/u", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
      
      # Python files and callback
      {"format": "PythonFile", "outputInterval": 1, "filename": "out/all/"+scenario_name, "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
      {"format": "PythonCallback", "outputInterval": 1, "filename": "out/all/"+scenario_name, "callback": handle_result_hyperelasticity, "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
    ],
    # 2. additional output writer that writes also the hydrostatic pressure
    "pressure": {   # output files for pressure function space (linear elements), contains pressure values, as well as displacements and velocities
      "OutputWriter" : [
        #{"format": "Paraview", "outputInterval": 1, "filename": "out/"+scenario_name+"/p", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
      ]
    },
    # 4. output writer for debugging, outputs files after each load increment, the geometry is not changed but u and v are written
    "LoadIncrements": {   
      "OutputWriter" : [
        #{"format": "Paraview", "outputInterval": 1, "filename": "out/load_increments", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
      ]
    },
  },
  "FiniteElementMethod" : {       # linear elasticity finite element method
    "meshName":             "3Dmesh_quadratic",           # mesh with quadratic Lagrange ansatz functions
    "inputMeshIsGlobal":    True,                         # boundary conditions are specified in global numberings, whereas the mesh is given in local numbering 
    "solverName":           "linearElasticitySolver",                   # reference to the linear solver
    "prefactor":            1.0,                                        # prefactor of the lhs, has no effect here
    "slotName":             "",
    "dirichletBoundaryConditions": elasticity_dirichlet_bc,             # the Dirichlet boundary conditions that define values for displacements u
    "dirichletOutputFilename":     None,                                # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
    "neumannBoundaryConditions":   elasticity_neumann_bc,               # Neumann boundary conditions that define traction forces on surfaces of elements
    "divideNeumannBoundaryConditionValuesByTotalArea": False,           # if the given Neumann boundary condition values under "neumannBoundaryConditions" are total forces instead of surface loads and therefore should be scaled by the surface area of all elements where Neumann BC are applied
    
    # material parameters
    "bulkModulus":          50,     # bulk modulus K, how much incompressible, high -> incompressible, low -> very compressible
    "shearModulus":         100,      # shear modulus, μ or G, "rigidity", how much shear stress response to shear deformation
    
    "OutputWriter" : [
      # Paraview files
      {"format": "Paraview", "outputInterval": 1, "filename": "out/"+scenario_name+"/u", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
      
      # Python files and callback
      {"format": "PythonFile", "outputInterval": 1, "filename": "out/all/"+scenario_name, "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
      {"format": "PythonCallback", "outputInterval": 1, "filename": "out/all/"+scenario_name, "callback": handle_result_linear_elasticity, "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
    ],
  },
  "NonlinearElasticitySolverFebio": {
    "durationLogKey": "febio",
    "tractionVector": traction_vector,                    # traction vector that is applied
    #"tractionElementNos": [(2*nz-1)*2*nx*2*ny + j*2*nx + i for j in range(2*ny) for i in range(2*nx)],    # elements on which traction is applied
    "tractionElementNos": [(nz-1)*nx*ny + j*nx + i for j in range(ny) for i in range(nx)],    # elements on which traction is applied
    "dirichletBoundaryConditionsMode": dirichlet_bc_mode, # "fix_all" or "fix_floating", how the bottom of the box will be fixed, fix_all fixes all nodes, fix_floating fixes all nodes only in z and the edges in x/y direction
    "materialParameters": material_parameters,            # c0, c1, k for Ψ = c0 * (I1-3) + c1 * (I2-3) + 1/2*k*(log(J))^2
    
    "meshName":             "3Dmesh_quadratic",           # mesh with quadratic Lagrange ansatz functions
    "inputMeshIsGlobal":    True,                         # boundary conditions are specified in global numberings, whereas the mesh is given in local numbering 
    "slotNames":            [],
    
    # 1. main output writer that writes output files using the quadratic elements function space. Writes displacements, velocities and PK2 stresses.
    "OutputWriter" : [
      
      # Paraview files
      {"format": "Paraview", "outputInterval": 1, "filename": "out/"+scenario_name+"/u", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
      
      # Python files and callback
      {"format": "PythonFile", "outputInterval": 1, "filename": "out/all/"+scenario_name, "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
      {"format": "PythonCallback", "outputInterval": 1, "filename": "out/all/"+scenario_name, "callback": handle_result_febio, "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
    ],
  },
}
