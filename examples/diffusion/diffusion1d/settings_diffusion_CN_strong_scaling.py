# Diffusion 1D

import sys

n = 1500
physical_extent = 15.0 # cm
dx = physical_extent/n # dx=0.01

dt = 0.001 # ms
end_time = 10 # 10000 iterations to encounter the cost of communication for long-term runs

Diffusion = 0.01 # almost comparable to shorten fast-twitch 
solver_type = "gmres"
preconditioner_type = "none"
scenario_name = "DiffusionImplicitCN"

#parse arguments  
print("{}".format(str(sys.argv)))

if len(sys.argv) >= 2:
  n = (int)(sys.argv[0]) # number of elements
print("n: {}".format(n))

if len(sys.argv) >= 3:
  scenario_name = sys.argv[1]
  print("scenario: [{}]".format(scenario_name))

if len(sys.argv) >= 4:
  solver_type = sys.argv[2]
  
if len(sys.argv) >= 5:
  preconditioner_type = sys.argv[3]

print("solver: [{}], preconditioner: [{}]".format(solver_type,preconditioner_type))

config = {
  "Meshes": {
    "mesh": {
      "nElements": n,
      "physicalExtent": physical_extent,
    }
  },
  "Solvers": {
    "implicitSolver": {
      "solverType": solver_type,
      "preconditionerType": preconditioner_type,
      "relativeTolerance": 1e-15,
    }
  },
  "scenarioName": scenario_name,
  "CrankNicolson" : {
     "timeStepWidth" : dt,
     "endTime": end_time,
     "initialValues": [20,20,40,50,20,20],
     "solverName": "implicitSolver",
     "durationLogKey": "duration_diffusion_1D_implicitSolver",    
     "meshName": "mesh",  
     "FiniteElementMethod" : {
        "diffusionTensor": Diffusion,
        "meshName": "mesh",
        "solverName": "implicitSolver",
     },
     "OutputWriter" : [
       #{"format": "Paraview", "outputInterval": 1, "filename": "out", "binaryOutput": "false", "fixedFormat": False},
       {"format": "PythonFile", "filename": "out/" + scenario_name, "outputInterval": 10000, "binary":False}
     ]
  }
}
