# isotropic Mooney Rivlin
import numpy as np
import sys, os

# number of elements
nx = 2    # 2
ny = 2    # 2
nz = 5    # 5

# boundary conditions (for quadratic elements)
dirichlet_bc = {}

xpos = 0.0
ypos = 0.0
zpos = 0.0

# fix bottom plane in z direction, displacements are quadratic
for j in range(0,2*ny+1):
  for i in range(0,2*nx+1):
    dirichlet_bc[j*(2*nx+1) + i] = [None,None,zpos,None,None,None]   # displacements and velocity

if False:
  # left plane
  for k in range(0,2*nz+1):
    for j in range(0,2*ny+1):
      dirichlet_bc[k*(2*nx+1)*(2*ny+1) + j*(2*nx+1)] = [xpos,None,None]

  # front plane
  for k in range(0,2*nz+1):
    for i in range(0,2*nx+1):
      dirichlet_bc[k*(2*nx+1)*(2*ny+1) + i] = [None,ypos,None]

  # vertical edge
  for k in range(0,2*nz+1):
    dirichlet_bc[k*(2*nx+1)*(2*ny+1)] = [xpos,ypos,None]

if True:
  # fix points on bottom horizontal edge in y and z direction
  for i in range(0,2*nx+1):
    dirichlet_bc[i] = [None,ypos,zpos,None,0]

if True:
  # horizontal edge
  for j in range(0,2*ny+1):
    dirichlet_bc[j*(2*nx+1)] = [xpos,None,zpos]

# fix corner completely
dirichlet_bc[0] = [xpos,ypos,zpos,None,0]

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
  for i in range(0,2*nx+1):
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

  print("factors = 1/36 * ")
  for j in range(0,2*ny+1):
    for i in range(0,2*nx+1):
      print("{}".format(36*factors[j * (2*nx+1) + i]),end =" ")
    print("")

  # compute force values
  total_force_x = 0
  total_force_y = 0
  total_force_z = 0
  
  # loop over nodes on top layer of geometry
  for j in range(0,2*ny+1):
    for i in range(0,2*nx+1):
      k = 2*nz+1-1
      factor = factors[j * (2*nx+1) + i]
      
      index = k*(2*nx+1)*(2*ny+1) + j*(2*nx+1) + i
      
      total_force_x += factor * s11_values[index]
      total_force_y += factor * s22_values[index]
      total_force_z += factor * s33_values[index]
  
  print("t: {}, force: {},{},{}".format(current_time, total_force_x, total_force_y, total_force_z))

# Neumann boundary conditions, specify upward force for top elements, slightly in y-direction
neumann_bc = [{"element": (nz-1)*nx*ny + j*nx + i, "constantVector": [0,1e-1,5e-1], "face": "2+"} for j in range(ny) for i in range(nx)]

dt = 1e-1
end_time = 100
output_interval = dt

#dirichlet_bc = {}
#neumann_bc = []

config = {
  "scenarioName": "3d_box",
  "DynamicHyperelasticitySolver": {
    #"numberTimeSteps": 1,
    "endTime": end_time,
    "timeStepWidth": dt,    
    "durationLogKey": "nonlinear",
    "timeStepOutputInterval": 1,
    
    #"materialParameters": [1.5,2.0],
    "materialParameters": [0.0,1.0],
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
    "physicalExtent": [nx, ny, nz],
    
    # nonlinear solver
    "relativeTolerance": 1e-10,         # 1e-10 relative tolerance of the linear solver
    "solverType": "preonly",            # type of the linear solver: cg groppcg pipecg pipecgrr cgne nash stcg gltr richardson chebyshev gmres tcqmr fcg pipefcg bcgs ibcgs fbcgs fbcgsr bcgsl cgs tfqmr cr pipecr lsqr preonly qcg bicg fgmres pipefgmres minres symmlq lgmres lcd gcr pipegcr pgmres dgmres tsirm cgls
    "preconditionerType": "lu",         # type of the preconditioner
    "maxIterations": 1e4,               # maximum number of iterations in the linear solver
    "snesMaxFunctionEvaluations": 1e8,  # maximum number of function iterations
    "snesMaxIterations": 50,            # maximum number of iterations in the nonlinear solver
    "snesRelativeTolerance": 1e-10,     # tolerance of the nonlinear solver
    
    #"dumpFilename": "out/r{}/m".format(sys.argv[-1]),   # dump system matrix and right hand side after every solve
    "dumpFilename": "",         # dump disabled
    "dumpFormat": "matlab",   # default, ascii, matlab
    
    # boundary and initial conditions
    "dirichletBoundaryConditions": dirichlet_bc,
    "neumannBoundaryConditions": neumann_bc,
    "updateDirichletBoundaryConditionsFunction": update_dirichlet_boundary_conditions,
    "updateDirichletBoundaryConditionsFunctionCallInterval": 1,
    
    "initialValuesDisplacements": [],
    #"initialValuesDisplacements": [[0.0,0.0,0.0] for i in range((2*nx+1)*(2*ny+1)*(2*nz+1-1))] + [[1.0,0.0,0.0] for i in range((2*nx+1)*(2*ny+1))],
    "initialValuesVelocities": [[0.01*z,0.0,0.0] for i in range((2*nx+1)*(2*ny+1)) for z in range((2*nz+1))],
    #"constantBodyForce": (1,0,0),     # e.g. for gravity
    
    "OutputWriter" : [   # output files for displacements function space (quadratic elements), contains displacements, velocities and PK2 stresses
      {"format": "Paraview", "outputInterval": 1, "filename": "out/u", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True},
      {"format": "PythonCallback", "outputInterval": 1, "callback": postprocess, "onlyNodalValues":True},
    ],
    "pressure": {   # output files for pressure function space (linear elements), contains pressure values, as well as displacements and velocities
      "OutputWriter" : [
        {"format": "Paraview", "outputInterval": 1, "filename": "out/p", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True},
      ]
    },
    "dynamic": {    # output of the dynamic solver, has additional virtual work values 
      "OutputWriter" : [   # output files for displacements function space (quadratic elements)
        #{"format": "Paraview", "outputInterval": int(output_interval/dt), "filename": "out/dynamic", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True},
        #{"format": "Paraview", "outputInterval": 1, "filename": "out/dynamic", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True},
      ],
    }
  }
}
