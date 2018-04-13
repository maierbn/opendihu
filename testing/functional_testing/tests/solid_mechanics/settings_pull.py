# 3D Incompressible Material
#
# command arguments: <name> <number elements>

import numpy as np
import scipy.integrate
import sys

name = ""

if len(sys.argv) > 0:
  if "check_results.py" not in sys.argv[0]:
    name = sys.argv[0]
    
    print("name: \"{}\"".format(name))

nx = 1
ny = 1
nz = 1
# [1,1,1] = 1 element, 8 dofs per element, 8 dofs, 24 unknowns

# dimensions: lx, lz
lx = 1.0
lz = 1.0
area_z = lz*lz
tmax = 2.0

dirichletBC = {
  0: 0.0, 1: 0.0, 2: 0.0,
  6: 0.0, 8: 0.0,
  12: 0.0, 13: 0.0,
  18: 0.0
} 

traction = {
  3: tmax/4.,
  9: tmax/4.,
  15: tmax/4.,
  21: tmax/4.
}

material_parameters = [6.352e-10, 3.627, 100]  # c0, c1, kappa

lambda = 4*material_parameters[1] / (tmax - 2*material_parameters[0])
analytical_pk2_stress = [2*material_parameters[0] + 4*material_parameters[1]/lambda, -2*material_parameters[1]/lambda, -2*material_parameters[1]/lambda]

print("lambda: ",lambda)

config = {
  "FiniteElementMethod" : {
    "nElements": [nx,ny,nz],
    "nodeDimension": 1,
    "physicalExtent": [lx,lz,lz],
    "DirichletBoundaryCondition": dirichletBC,  # displacement Dirichlet bc
    "relativeTolerance": 1e-15,
    "rightHandSide": traction,  # surface traction or body force
    "materialParameters": material_parameters,  # c0, c1, kappa
    "analyticJacobian": False,   # False = compute Jacobian by finite differences
  },
  "OutputWriter" : [
    #{"format": "Paraview", "outputInterval": 1, "filename": "out", "binaryOutput": "false", "fixedFormat": False},
    #{"format": "ExFile", "filename": "out/"+name, "outputInterval": 2},
    {"format": "PythonFile", "filename": "out/"+name, "outputInterval": 5, "binary":False, "onlyNodalValues":True},
  ]
}

# output config in a readable format
if True:
  import pprint 
  pp = pprint.PrettyPrinter()
  pp.pprint(config)
