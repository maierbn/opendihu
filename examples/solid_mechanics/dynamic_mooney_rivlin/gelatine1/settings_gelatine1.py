# Dynamic, isotropic Mooney Rivlin, mixed formulation (Taylor-Hood)
#
# 
import numpy as np
import sys, os

# Define number of elements
# The number of nodes in (2*nx+1) x (2*ny+1) x (2*nz+1) because we have a quadratic mesh
# z axis is pointing upwards
# y axis is pointing in moving direction
nx = 4
ny = 2
nz = 6

# set up Dirichlet boundary condition values (indexing is for quadratic elements)
dirichlet_bc = {}
xpos = 0.0
ypos = 0.0
zpos = 0.0

# Fix top plane, set y velocity to 0
for j in range(0,2*ny+1):
  for i in range(0,2*nx+1):
    k = 2*nz
    dirichlet_bc[k*(2*nx+1)*(2*ny+1) + j*(2*nx+1) + i] = [xpos, ypos, zpos, None, 0.0, None]   # displacements [ux,uy,uz] and velocity [vx,vy,vz]

# Function to update dirichlet boundary conditions over time
# This function returns "dirichlet_bc". Only those entries can be updated that were also initially set.
def update_dirichlet_boundary_conditions(t):
  
  # parameters for movement
  duration = 5.0
  amplitude = 1.0
  
  # define movement in y-direction
  def y(t):
    t *= np.pi/duration
    if t > np.pi:
      return 2*amplitude
    return (1 + np.cos(t - np.pi))*amplitude
  
  # define velocity of movement, dy(t)/dt
  def dydt(t):
    t *= np.pi/duration
    if t > np.pi:
      return 0
    return -np.sin(t - np.pi)*amplitude
    
  # set d
  for j in range(0,2*ny+1):
    for i in range(0,2*nx+1):
      k = 2*nz
      dirichlet_bc[k*(2*nx+1)*(2*ny+1) + j*(2*nx+1) + i] = [xpos, ypos+y(t), zpos, None, dydt(t), None]  # displacements [ux,uy,uz] and velocity [vx,vy,vz]

  return dirichlet_bc

# Function to postprocess the output
# This function gets periodically called by the running simulation. 
# It provides all current variables for each node: geometry (position), u, v, stress
def postprocess(result):
  result = result[0]
  # print result for debugging
  #print(result)
  
  # get current time
  current_time = result["currentTime"]
  timestep_no = result["timeStepNo"]
  
  # parse variables
  field_variables = result["data"]
  
  # output all field variables for debugging
  for i,f in enumerate(field_variables):
    print(i,field_variables[i]["name"])
  
  # field_variables[0] is the geometry
  # field_variables[1] is the displacements u
  # field_variables[2] is the velocities v
  # field_variables[3] is the PK2-Stress (Voigt)
  # field_variables[4] is the material traction
  
  traction_components = field_variables[3]["components"]
  
  # traction values contains the traction vector in reference configuration
  traction1_values = traction_components[0]["values"]   # traction in x-direction
  traction2_values = traction_components[1]["values"]   # traction in y-direction
  traction3_values = traction_components[2]["values"]   # traction in z-direction

  # Integrate total forces, e.g. Fx = ∫∫ S_11(x,y) • ϕ(x,y) dxdy. The integration domain is the top layer of the geometry.
  # The result corresponds to the bearing forces.
  # Integration is done by adding the values with factors

  # Integration stencils for ∫∫ ϕ_i dxdy with ϕ quadratic Lagrange function
  # nodes i of a single quadratic 2D element:
  #   1/36  1/9  1/36      primary   secondary  primary
  #   1/9   4/9  1/9   =   secondary tertiary   secondary
  #   1/36  1/9  1/36      primary   secondary  primary
  #
  # values in total mesh: 
  # interior primary node: 4*1/36 = 1/9, interior secondary node: 2*1/9 = 2/9, interior tertiary node: 4/9
  # edge primary node: 2*1/36 = 1/18, edge secondary node: 1/9
  # corner primary node: 1/36

  # store all factors
  factors = [1./9. for _ in range((2*nx+1) * (2*ny+1))]
  
  # set factors for secondary nodes
  for j in range(2*ny+1):
    for i in range(2*nx+1):
      if (i+j) % 2 == 1:
        factors[j * (2*nx+1) + i] = 2./9.
      
  # set factors for tertiary nodes
  for j in range(1,2*ny+1,2):
    for i in range(1,2*nx+1,2):
      factors[j * (2*nx+1) + i] = 4./9.

  # edges, primary nodes
  for i in range(0,2*nx+1,2):
    factors[(2*ny) * (2*nx+1) + i] = 1./18.
    factors[i] = 1./18.
  for j in range(0,2*ny+1,2):
    factors[j * (2*nx+1) + 0] = 1./18.
    factors[j * (2*nx+1) + (2*nx)] = 1./18.
    
  # edges, secondary nodes
  for i in range(1,2*nx+1,2):
    factors[(2*ny) * (2*nx+1) + i] = 1./9.
    factors[i] = 1./9.
  for j in range(1,2*ny+1,2):
    factors[j * (2*nx+1) + 0] = 1./9.
    factors[j * (2*nx+1) + (2*nx)] = 1./9.
    
  # corners
  factors[0] = 1./36
  factors[2*nx] = 1./36
  factors[(2*ny) * (2*nx+1) + 0] = 1./36
  factors[(2*ny) * (2*nx+1) + (2*nx)] = 1./36

  # output factors for debugging
  if True:
    print("factors = 1/36 * ")
    for j in range(0,2*ny+1):
      for i in range(0,2*nx+1):
        print("{}".format(36*factors[j * (2*nx+1) + i]),end =" ")
      print("")

  # integrate force values
  total_force_x = 0
  total_force_y = 0
  total_force_z = 0
  
  # loop over nodes on top layer of geometry
  for j in range(0,2*ny+1):
    for i in range(0,2*nx+1):
      k = 2*nz+1-1
      factor = factors[j * (2*nx+1) + i]
      
      index = k*(2*nx+1)*(2*ny+1) + j*(2*nx+1) + i
      
      total_force_x += factor * traction1_values[index]
      total_force_y += factor * traction2_values[index]
      total_force_z += factor * traction3_values[index]
  
  # output values
  print("t: {}, force: {},{},{}".format(current_time, total_force_x, total_force_y, total_force_z))

dt = 0.5
end_time = 100
output_interval = dt

neumann_bc = []

config = {
  "scenarioName": "3d_box",               # specifier to find simulation run in log file
  "logFormat":    "csv",                  # "csv" or "json", format of the lines in the log file, csv gives smaller files
  "solverStructureDiagramFile":     "solver_structure.txt",     # output file of a diagram that shows data connection between solvers
  "mappingsBetweenMeshesLogFile": "mappings_between_meshes_log.txt",    # log file for mappings 
  "DynamicHyperelasticitySolver": {
    #"numberTimeSteps": 1,
    "endTime": end_time,                  # end time of simulation
    "timeStepWidth": dt,                  # time step width 
    "durationLogKey": "nonlinear",        # key to find duration of this solver in the log file
    "timeStepOutputInterval": 1,          # how often the current time step should be printed to console
    
    #"materialParameters": [1.5,2.0],
    "materialParameters": [0.0,1.0],      # material parameters of the Mooney-Rivlin material
    "density": 1.0,                       # density of the material
    "displacementsScalingFactor": 1.0,    # scaling factor for displacements, only set to sth. other than 1 only to increase visual appearance for very small displacements
    "residualNormLogFilename": "log_residual_norm.txt",   # log file where residual norm values of the nonlinear solver will be written
    "useAnalyticJacobian": True,          # whether to use the analytically computed jacobian matrix in the nonlinear solver (fast)
    "useNumericJacobian": False,          # whether to use the numerically computed jacobian matrix in the nonlinear solver (slow), only works with non-nested matrices, if both numeric and analytic are enable, it uses the analytic for the preconditioner and the numeric as normal jacobian
      
    "dumpDenseMatlabVariables": False,    # whether to have extra output of matlab vectors, x,r, jacobian matrix (very slow)
    # if useAnalyticJacobian,useNumericJacobian and dumpDenseMatlabVariables all all three true, the analytic and numeric jacobian matrices will get compared to see if there are programming errors for the analytic jacobian
    
    # mesh
    "nElements": [nx, ny, nz],            # number of elements in x,y and z direction
    "inputMeshIsGlobal": True,            # if all indices are global (set to True)
    "physicalExtent": [nx, ny, nz],       # physical size in [m]
    "physicalOffset": [0, 0, 0],          # bottom front left position of the mesh
    
    # nonlinear solver
    "relativeTolerance":          1e-10,                        # 1e-10 relative tolerance of the linear solver
    "absoluteTolerance":          1e-10,                        # 1e-10 absolute tolerance of the residual of the linear solver       
    "solverType":                 "preonly",                    # type of the linear solver: cg groppcg pipecg pipecgrr cgne nash stcg gltr richardson chebyshev gmres tcqmr fcg pipefcg bcgs ibcgs fbcgs fbcgsr bcgsl cgs tfqmr cr pipecr lsqr preonly qcg bicg fgmres pipefgmres minres symmlq lgmres lcd gcr pipegcr pgmres dgmres tsirm cgls
    "preconditionerType":         "lu",                         # type of the preconditioner
    "maxIterations":              1e4,                          # maximum number of iterations in the linear solver
    "snesMaxFunctionEvaluations": 1e8,                          # maximum number of function iterations
    "snesMaxIterations":          10,                           # maximum number of iterations in the nonlinear solver
    "snesRelativeTolerance":      1e-5,                         # relative tolerance of the nonlinear solver
    "snesLineSearchType":         "l2",                         # type of linesearch, possible values: "bt" "nleqerr" "basic" "l2" "cp" "ncglinear"
    "snesAbsoluteTolerance":      1e-5,                         # absolute tolerance of the nonlinear solver
    "snesRebuildJacobianFrequency": 10,                         # how often the jacobian should be recomputed, -1 indicates NEVER rebuild, 1 means rebuild every time the Jacobian is computed within a single nonlinear solve, 2 means every second time the Jacobian is built etc. -2 means rebuild at next chance but then never again 
    
    "dumpFilename": "",                                         # filename of vector and matrix dump, "" means disabled
    "dumpFormat": "matlab",                                     # format of vector and matrix dump
    
    #"loadFactors":                [0.5, 1.0],                   # load factors for every timestep
    "loadFactors":                [],                           # no load factors, solve problem directly
    "loadFactorGiveUpThreshold": 0.1,                           # if the adaptive time stepping produces a load factor smaller than this value, the solution will be accepted for the current timestep, even if it did not converge fully to the tolerance
    "nNonlinearSolveCalls":       1,                            # how often the nonlinear solve should be called
    
    # boundary and initial conditions
    "dirichletBoundaryConditions": dirichlet_bc,                        # the initial Dirichlet boundary conditions that define values for displacements u and velocity v
    "neumannBoundaryConditions":   neumann_bc,                          # Neumann boundary conditions that define traction forces on surfaces of elements
    "divideNeumannBoundaryConditionValuesByTotalArea": True,            # if the given Neumann boundary condition values under "neumannBoundaryConditions" are total forces instead of surface loads and therefore should be scaled by the surface area of all elements where Neumann BC are applied
    
    "updateDirichletBoundaryConditionsFunction":  update_dirichlet_boundary_conditions,  # function that updates the Dirichlet BCs while the simulation is running
    "updateDirichletBoundaryConditionsFunctionCallInterval": 1,                          # every which step the update function should be called, 1 means every time step
    "updateNeumannBoundaryConditionsFunction":    None,                                  # function that updates the Neumann BCs while the simulation is running
    "updateNeumannBoundaryConditionsFunctionCallInterval": 1,                            # every which step the update function should be called, 1 means every time step
    
    "initialValuesDisplacements": [[0.0,0.0,0.0] for _ in range((2*nx+1) * (2*ny+1) * (2*nz+1))],     # the initial values for the displacements, vector of values for every node [[node1-x,y,z], [node2-x,y,z], ...]
    "initialValuesVelocities":    [[0.0,0.0,0.0] for _ in range((2*nx+1) * (2*ny+1) * (2*nz+1))],     # the initial values for the velocities, vector of values for every node [[node1-x,y,z], [node2-x,y,z], ...]
    "extrapolateInitialGuess":     True,                                # if the initial values for the dynamic nonlinear problem should be computed by extrapolating the previous displacements and velocities
    "constantBodyForce": (0,0,-1e-1),                 # a constant force that acts on the whole body, e.g. for gravity
    
    "dirichletOutputFilename":     "out/dirichlet_boundary_conditions",                                # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
    
    # define which file formats should be written
    # 1. main output writer that writes output files using the quadratic elements function space. Writes displacements, velocities and PK2 stresses.
    "OutputWriter" : [
      
      # Paraview files
      {"format": "Paraview", "outputInterval": 1, "filename": "out/u", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
      
      # Python callback function "postprocess" (at the top of this scritpt)
      {"format": "PythonCallback", "outputInterval": 1, "callback": postprocess, "onlyNodalValues":True, "filename": "", "fileNumbering": "incremental"},
    ],
    # 2. additional output writer that writes also the hydrostatic pressure
    "pressure": {   # output files for pressure function space (linear elements), contains pressure values, as well as displacements and velocities
      "OutputWriter" : [
        #{"format": "Paraview", "outputInterval": 1, "filename": "out/p", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
      ]
    },
    # 3. additional output writer that writes virtual work terms
    "dynamic": {    # output of the dynamic solver, has additional virtual work values 
      "OutputWriter" : [   # output files for displacements function space (quadratic elements)
        #{"format": "Paraview", "outputInterval": int(output_interval/dt), "filename": "out/dynamic", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
        #{"format": "Paraview", "outputInterval": 1, "filename": "out/dynamic", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
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
