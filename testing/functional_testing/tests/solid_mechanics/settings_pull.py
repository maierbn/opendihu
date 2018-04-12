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
    
    print "name: \"{}\"".format(name)

nx = 1
ny = 1
nz = 1
# [1,1,1] = 1 element, 8 dofs per element, 8 dofs, 24 unknowns

dirichletBC = {
  0: 0.0, 1: 0.0, 2: 0.0,
  6: 0.0, 8: 0.0,
  12: 0.0, 13: 0.0,
  18: 0.0
} 

traction = {
  3: 1.0,
  9: 1.0,
  15: 1.0,
  21: 1.0
}

config = {
  "FiniteElementMethod" : {
    "nElements": [nx,ny,nz],
    "nodeDimension": 1,
    "physicalExtent": [2.*nx,2.*ny,2.*nz],
    "DirichletBoundaryCondition": dirichletBC,  # displacement Dirichlet bc
    "relativeTolerance": 1e-15,
    "rightHandSide": traction,  # surface traction or body force
    "materialParameters": [6.352e-10, 3.627, 100],  # c0, c1, kappa
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
