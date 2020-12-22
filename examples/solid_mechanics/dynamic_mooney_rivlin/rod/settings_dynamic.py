# isotropic Mooney Rivlin
import numpy as np
import sys, os

# number of elements
nx = 2    # 2
ny = 2    # 2
nz = 3    # 5

# number of nodes
mx = 2*nx + 1
my = 2*ny + 1
mz = 2*nz + 1

# boundary conditions (for quadratic elements)
dirichlet_bc = {}

xpos = 0.0
ypos = 0.0
zpos = 0.0

# fix bottom plane in z direction, displacements are quadratic
for j in range(0,my):
  for i in range(0,mx):
    dirichlet_bc[j*mx + i] = [None,None,zpos,None,None,None]   # displacements and velocity

if False:
  # left plane
  for k in range(0,mz):
    for j in range(0,my):
      dirichlet_bc[k*mx*my + j*mx] = [xpos,None,None]

  # front plane
  for k in range(0,mz):
    for i in range(0,mx):
      dirichlet_bc[k*mx*my + i] = [None,ypos,None]

  # vertical edge
  for k in range(0,mz):
    dirichlet_bc[k*mx*my] = [xpos,ypos,None]

if True:
  # fix points on bottom horizontal edge in y and z direction
  for i in range(0,mx):
    dirichlet_bc[i] = [None,ypos,zpos,None,0,None]

if True:
  # horizontal edge
  for j in range(0,my):
    dirichlet_bc[j*mx] = [xpos,None,zpos,None,None,None]

# fix corner completely
dirichlet_bc[0] = [xpos,ypos,zpos,None,0,None]

# function used to update dirichlet boundary conditions
def update_dirichlet_boundary_conditions(t):
  
  def y(t):
    a = 1.0
    t *= 0.1
    if t > np.pi:
      return 2*a
    return (1 + np.cos(t - np.pi))*a
  
  def dydt(t):
    a = 1.0
    t *= 0.1
    if t > np.pi:
      return 0
    return (-np.sin(t - np.pi))*a
  
  # fix points on bottom horizontal edge in y and z direction
  for i in range(0,mx):
    dirichlet_bc[i] = [None,ypos+y(t),zpos,None,dydt(t)]
    
  # fix corner completely
  dirichlet_bc[0] = [xpos,ypos+y(t),zpos,None,dydt(t)]
  
  return dirichlet_bc

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
  stress_components = field_variables[3]["components"]
  # stress_components contains the symmetric 3x3 material stress tensor, in this order: S_11, S_22, S_33, S_12, S_13, S_23

  s11_values = stress_components[0]["values"]
  s22_values = stress_components[1]["values"]
  s33_values = stress_components[2]["values"]

  mx = 2*result["nElementsLocal"][0] + (1 if result["hasFullNumberOfNodes"][0] else 0)
  my = 2*result["nElementsLocal"][1] + (1 if result["hasFullNumberOfNodes"][1] else 0)
  mz = 2*result["nElementsLocal"][2] + (1 if result["hasFullNumberOfNodes"][2] else 0)

  # integrate total forces

  # integration stencils for ∫∫ ϕ_i dxdy with ϕ quadratic Lagrange function
  # nodes i of a single quadratic 2D element:
  # 1/36  1/9  1/36      primary   secondary  primary
  # 1/9   4/9  1/9   =   secondary tertiary   secondary
  # 1/36  1/9  1/36      primary   secondary  primary
  #
  # total mesh: 
  # interior primary node: 4*1/36 = 1/9, interior secondary node: 2*1/9 = 2/9, interior tertiary node: 4/9
  # edge primary node: 2*1/36 = 1/18, edge secondary node: 1/9
  # corner primary node: 1/36

  factors = [1./9. for _ in range(mx * my)]
  
  # set factors for secondary nodes
  for j in range(my):
    for i in range(mx):
      if (i+j) % 2 == 1:
        factors[j * mx + i] = 2./9.
      
  # set factors for tertiary nodes
  for j in range(1,my,2):
    for i in range(1,mx,2):
      factors[j * mx + i] = 4./9.

  # edges, primary nodes
  for i in range(0,mx,2):
    factors[(my-1) * mx + i] = 1./18.
    factors[i] = 1./18.
  for j in range(0,my,2):
    factors[j * mx + 0] = 1./18.
    factors[j * mx + (mx-1)] = 1./18.
    
  # edges, secondary nodes
  for i in range(1,mx,2):
    factors[(my-1) * mx + i] = 1./9.
    factors[i] = 1./9.
  for j in range(1,my,2):
    factors[j * mx + 0] = 1./9.
    factors[j * mx + (mx-1)] = 1./9.
    
  # corners
  factors[0] = 1./36
  factors[2*nx] = 1./36
  factors[(my-1) * mx + 0] = 1./36
  factors[(my-1) * mx + (mx-1)] = 1./36

  if False:
    print("factors = 1/36 * ")
    for j in range(0,my):
      for i in range(0,mx):
        print("{}".format(36*factors[j * mx + i]),end =" ")
      print("")

  # compute force values
  total_force_x = 0
  total_force_y = 0
  total_force_z = 0
  
  # loop over nodes on top layer of geometry
  for j in range(0,my):
    for i in range(0,mx):
      k = mz-1
      factor = factors[j * mx + i]
      
      index = k*mx*my + j*mx + i
      total_force_x += factor * s11_values[index]
      total_force_y += factor * s22_values[index]
      total_force_z += factor * s33_values[index]
  
  print("t: {}, total force on top bearing: {},{},{}".format(current_time, total_force_x, total_force_y, total_force_z))

# Neumann boundary conditions, specify upward force for top elements, slightly in y-direction
neumann_bc = [{"element": (nz-1)*nx*ny + j*nx + i, "constantVector": [0,1e-1,1.0], "face": "2+"} for j in range(ny) for i in range(nx)]

# fiber directions
fiber_meshes = {}
fiber_mesh_names = []

# each fiber will have mz nodes
n_elements_fiber = mz-1

for j in range(my):
  for i in range(mx):
    fiber_no = j*nx + i
    
    # determine start position of fiber in (x,y)-plane
    x = 0.5 + i
    y = -2.0 + 0.2 + j
    angle = 20./180.*np.pi
    
    # loop over points of a single fiber
    node_positions = []
    for z in range(mz):
        
      h = 3.0*z/nz
      x_pos = x
      y_pos = y + np.sin(angle)*h
      z_pos = 0.0 + np.cos(angle)*h
      node_positions.append([x_pos,y_pos,z_pos])
    
    mesh_name = "fiber{}".format(fiber_no)
    fiber_mesh_names.append(mesh_name)
    
    fiber_meshes[mesh_name] = {
      "nodePositions": node_positions,
      "nElements": [n_elements_fiber],
      "inputMeshIsGlobal": True,
      "nRanks": [1],
    }
    
# load
constant_body_force = (0,0,0)

# time parameters
# -----------------
dt = 1e-1
end_time = 100
output_interval = dt

#dirichlet_bc = {}
#neumann_bc = []

config = {
  "scenarioName": "dynamic_rod",
  "logFormat":    "csv",     # "csv" or "json", format of the lines in the log file, csv gives smaller files
  "solverStructureDiagramFile":     "solver_structure.txt",     # output file of a diagram that shows data connection between solvers
  "mappingsBetweenMeshesLogFile":   "mappings_between_meshes.txt",   # log file for mappings between meshes
  
  "Meshes": fiber_meshes,  
  "Solvers": {
    "nonlinearSolver": {
      # nonlinear solver
      "relativeTolerance": 1e-5,         # 1e-10 relative tolerance of the linear solver
      "absoluteTolerance": 1e-5,         # 1e-10 absolute tolerance of the residual of the linear solver    
      "solverType": "preonly",            # type of the linear solver: cg groppcg pipecg pipecgrr cgne nash stcg gltr richardson chebyshev gmres tcqmr fcg pipefcg bcgs ibcgs fbcgs fbcgsr bcgsl cgs tfqmr cr pipecr lsqr preonly qcg bicg fgmres pipefgmres minres symmlq lgmres lcd gcr pipegcr pgmres dgmres tsirm cgls
      "preconditionerType": "lu",         # type of the preconditioner
      "maxIterations": 1e4,               # maximum number of iterations in the linear solver
      "dumpFilename": "",#"out/m",            # filename for output of solver matrix
      "dumpFormat": "matlab",             # default, ascii, matlab
      "snesMaxFunctionEvaluations": 1e8,  # maximum number of function iterations
      "snesMaxIterations": 50,             # maximum number of iterations in the nonlinear solver
      "snesRebuildJacobianFrequency": 5,  # frequency with which the jacobian is newly computed
      "snesRelativeTolerance": 1e-5,      # relative tolerance of the nonlinear solver
      "snesLineSearchType": "l2",         # type of linesearch, possible values: "bt" "nleqerr" "basic" "l2" "cp" "ncglinear"
      "snesAbsoluteTolerance": 1e-5,      # absolute tolerance of the nonlinear solver
    }
  },
  "DynamicHyperelasticitySolver": {
    #"numberTimeSteps": 1,
    "endTime": end_time,
    "timeStepWidth": dt,    
    "durationLogKey": "nonlinear",
    "timeStepOutputInterval": 1,
    
    #"materialParameters": [1.5,2.0],
    "materialParameters": [0.0,1.0,2.0,3.0],
    "density": 1.0,
    "displacementsScalingFactor": 1.0,   # scaling factor for displacements, set to sth. other than 1 only to increase visual appearance for very small displacements
    "residualNormLogFilename": "log_residual_norm.txt",
    "useAnalyticJacobian": True,
    "useNumericJacobian": False,   # only works with non-nested matrices, if both numeric and analytic are enable, it uses the analytic for the preconditioner and the numeric as normal jacobian
      
    "dumpDenseMatlabVariables": False,   # extra output of matlab vectors, x,r, jacobian matrix, very slow
    # if useAnalyticJacobian,useNumericJacobian and dumpDenseMatlabVariables all all three true, the analytic and numeric jacobian matrices will get compared to see if there are programming errors for the analytic jacobian
    
    # mesh
    "nElements": [nx, ny, nz],
    "inputMeshIsGlobal": True,
    "physicalExtent": [2.0, 2.0, 3.0],
    "physicalOffset": [0, 0, 0],        # offset/translation where the whole mesh begins
    
    "fiberMeshNames": fiber_mesh_names,   # fiber meshes that will be used to determine the fiber direction
    
    # nonlinear solver
    "solverName": "nonlinearSolver",
    
    #"dumpFilename": "out/r{}/m".format(sys.argv[-1]),   # dump system matrix and right hand side after every solve
    "dumpFilename": "",         # dump disabled
    "dumpFormat": "matlab",   # default, ascii, matlab
    
    #"loadFactors":  [0.1, 0.2, 0.35, 0.5, 1.0],   # load factors for every timestep
    "loadFactors": [],                 # no load factors, solve problem directly
    "nNonlinearSolveCalls": 1,         # how often the nonlinear solve should be repeated
    "loadFactorGiveUpThreshold": 0.1,   # if the adaptive time stepping produces a load factor smaller than this value, the solution will be accepted for the current timestep, even if it did not converge fully to the tolerance
    
    # boundary and initial conditions
    "dirichletBoundaryConditions": dirichlet_bc,
    "neumannBoundaryConditions": neumann_bc,
    "divideNeumannBoundaryConditionValuesByTotalArea": True,            # if the given Neumann boundary condition values under "neumannBoundaryConditions" are total forces instead of surface loads and therefore should be scaled by the surface area of all elements where Neumann BC are applied
    #"updateDirichletBoundaryConditionsFunction": update_dirichlet_boundary_conditions,
    "updateDirichletBoundaryConditionsFunction": None,
    "updateDirichletBoundaryConditionsFunctionCallInterval": 1,
    
    "initialValuesDisplacements": [[0.0,0.0,0.0] for i in range(mx*my*mz)],
    #"initialValuesDisplacements": [[0.0,0.0,0.0] for i in range(mx*my*(mz-1))] + [[1.0,0.0,0.0] for i in range(mx*my)],
    "initialValuesVelocities": [[0.1*z,0.0,0.0] for z in range(mz) for i in range(mx*my)],
    "extrapolateInitialGuess":    True,                     # if the initial values for the dynamic nonlinear problem should be computed by extrapolating the previous displacements and velocities
    "constantBodyForce": constant_body_force,     # e.g. for gravity
    
    "dirichletOutputFilename":     "out/dirichlet_boundary_conditions",                                # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
    
    "OutputWriter" : [   # output files for displacements function space (quadratic elements), contains displacements, velocities and PK2 stresses
      {"format": "Paraview", "outputInterval": 5, "filename": "out/u", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
      {"format": "PythonCallback", "outputInterval": 5, "callback": postprocess, "onlyNodalValues":True, "filename": "", "fileNumbering": "incremental"},
    ],
    "pressure": {   # output files for pressure function space (linear elements), contains pressure values, as well as displacements and velocities
      "OutputWriter" : [
        {"format": "Paraview", "outputInterval": 5, "filename": "out/p", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
      ]
    },
    "dynamic": {    # output of the dynamic solver, has additional virtual work values 
      "OutputWriter" : [   # output files for displacements function space (quadratic elements)
        #{"format": "Paraview", "outputInterval": int(output_interval/dt), "filename": "out/dynamic", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
        #{"format": "Paraview", "outputInterval": 1, "filename": "out/dynamic", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
      ],
    },
    # output writer for debugging, outputs files after each load increment, the geometry is not changed but u and v are written
    "LoadIncrements": {   
      "OutputWriter" : [
        {"format": "Paraview", "outputInterval": 5, "filename": "out_static/p", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
      ]
    },
  }
}
