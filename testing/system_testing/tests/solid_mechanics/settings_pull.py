# 3D Incompressible Material
#
# command arguments: <name> <number elements>

# reset && sdd && ./mooney_rivlin_incompressible_penalty ../settings_pull.py -snes_type newtonls -snes_linesearch_type basic  -snes_linesearch_damping 0.0001 -snes_max_it 100000  -snes_max_funcs 100000

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
tmax = 0.5
that = 0.5

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

#material_parameters = [6.352e-10, 3.627, 0]  # c0, c1, kappa
material_parameters = [1.0, 0.0, 100.0]  # c0, c1, kappa
# kappa = 0 means disable incompressibility and volume effects


#print("lambdaValue: ",lambdaValue,"> 1")

config = {
  "FiniteElementMethod" : {
    "nElements": [nx,ny,nz],
    "nodeDimension": 1,
    "physicalExtent": [lx,lz,lz],
    "dirichletBoundaryCondition": dirichletBC,  # displacement Dirichlet bc
    #"tractionReferenceConfiguration": [
    #  {"element": 0, "face": "0+", "constantValue": tmax/(lz**2)}   # face can be one of "0+", "0-", "1+", "1-", "2+", "2-". dofValues are the element-local node/dof numbers (e.g. 0-26 for 3D quadratic Lagrange elements), the entries are vectors of dimension D with respect to the local element coordinate system
    #],
    "tractionCurrentConfiguration": [
      {"element": 0, "face": "0+", "constantValue": that}   # face can be one of "0+", "0-", "1+", "1-", "2+", "2-". dofValues are the element-local node/dof numbers (e.g. 0-26 for 3D quadratic Lagrange elements), the entries are vectors of dimension D with respect to the local element coordinate system
    ],
    #"tractionReferenceConfiguration": [
    #  {"element": 1, "face": "0+", "dofVectors": {0: [tmax,0,0], 1: [tmax,0,0], 2: [tmax,0,0], 3: [tmax,0,0]}},  # face can be one of "0+", "0-", "1+", "1-", "2+", "2-". dofValues are the element-local node/dof numbers (e.g. 0-26 for 3D quadratic Lagrange elements), the entries are vectors of dimension D with respect to the local element coordinate system
    #  {"element": 1, "face": "1-", "dofVectors": {5: [1,2,3], 7:[8,3,4]}}, 
    #  {"element": 3, "face": "0-", "constantVector": tmax},   # a constant vector for constant traction on that element
    #  {"element": 3, "face": "0+", "constantValue": 5.12},   # a single constant value for constant traction on that element, in direction of outward normal
    #],
    #"tractionCurrentConfiguration": [
    #  
    #],
    #"bodyForceReferenceConfiguration": {0: [tmax,0,0], 5: [tmax,tmax,tmax]},   # {<element global no.>: <force vector, constant in element>, ...}
    #"bodyForceCurrentConfiguration": {},
    "relativeTolerance": 1e-15,
    "rightHandSide": traction,  # surface traction T or body force B, both in material description
    "materialParameters": material_parameters,  # c0, c1, kappa
    "analyticJacobian": False,   # False = compute Jacobian by finite differences
  },
  "OutputWriter" : [
    #{"format": "Paraview", "outputInterval": 1, "filename": "out", "binary": "false", "fixedFormat": False},
    #{"format": "ExFile", "filename": "out/"+name, "outputInterval": 2},
    {"format": "PythonFile", "filename": "out/"+name, "outputInterval": 5, "binary":False, "onlyNodalValues":True},
  ]
}

# output config in a readable format
if True:
  import pprint 
  pp = pprint.PrettyPrinter()
  pp.pprint(config)
