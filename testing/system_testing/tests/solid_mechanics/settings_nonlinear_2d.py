import random

nx = 10
ny = 5

dirichletBC = {
  0: 0.0, 1: 0.0,
} 

for y in range(1+2*ny):
  dirichletBC[y*(2 + 4*nx)] = 0.0

material_parameters = [5.0, 3.0, 100.0]  # c0, c1, kappa

random.seed(0)
node_positions = []
for y in range(1+2*ny):
  for x in range(1+2*nx):
    node_positions.append([x+random.random()*0.5,y+random.random()*0.5])

traction = []
for y in range(ny):
  element_no = nx-1 + y*nx
  traction.append({"element": element_no, "face": "0+", "constantValue": 0.1})
  print("traction on el ",element_no,",y:",y)


config = {
  "logFormat": "csv",
  "FiniteElementMethod" : {
    "nElements": [nx,ny],
    "physicalExtent": [10.0,10.0,10.0],
    "nodePositions": node_positions,
    "dirichletBoundaryCondition": dirichletBC,  # displacement Dirichlet bc
    "tractionReferenceConfiguration": traction,
    "relativeTolerance": 1e-15,
    "materialParameters": material_parameters,  # c0, c1, kappa
    "analyticJacobian": True,
    "numericJacobian": True,
    "logfile": "residual_norm.txt",
      
    "OutputWriter" : [
      {"format": "PythonFile", "filename": "out/nonlinear_test", "outputInterval": 1, "binary":False, "onlyNodalValues":True},
    ],
    "outputIntermediateSteps": True
  },
}
