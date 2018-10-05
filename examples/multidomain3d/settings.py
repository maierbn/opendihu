import numpy as np
import pickle

# global parameters
PMax = 7.3              # maximum stress [N/cm^2]
Conductivity = 3.828    # sigma, conductivity [mS/cm]
Am = 500.0              # surface area to volume ratio [cm^-1]
Cm = 0.58               # membrane capacitance [uF/cm^2]
innervation_zone_width = 1.  # cm
innervation_zone_width = 0.  # cm
mesh_file = "../input/mesh"

# load mesh
with open(mesh_file, "rb") as f:
  mesh_data = pickle.load(f, encoding='latin1')
#
#  "node_positions": node_positions, 
#  "linear_elements": linear_elements, 
#  "quadratic_elements": quadratic_elements, 
#  "seed_points": seed_points,
#  "bottom_nodes": bottom_node_indices,
#  "top_nodes": top_node_indices,
#  "n_linear_elements_per_coordinate_direction": n_linear_elements_per_coordinate_direction,
#  "n_quadratic_elements_per_coordinate_direction": n_quadratic_elements_per_coordinate_direction,
#


# boundary conditions
potential_flow_bc = {}
for bottom_node_index in mesh_data["bottom_nodes"]:
  potential_flow_bc[bottom_node_index] = 0.0
for top_node_index in mesh_data["top_nodes"]:
  potential_flow_bc[top_node_index] = 1.0
  
config = {
  "Meshes": {
    "mesh": {
      "nElements": mesh_data["n_linear_elements_per_coordinate_direction"],
      "nodePositions": mesh_data["node_positions"],
      "inputMeshIsGlobal": True,
      "setHermiteDerivatives": False
    }
  },
  "Solvers": {
    "potentialFlowSolver": {
      "relativeTolerance": 1e-15,
    },
    "activationSolver": {
      "relativeTolerance": 1e-15,
    }
  },
  "MultidomainSolver" : {
    "nCompartments": 1,
    "am": Am,
    "cm": Cm,
    "timeStepWidth": 1e-3,
    "endTime": 1.0,
    
    "PotentialFlow": {
      "FiniteElementMethod" : {  
        "meshName": "mesh",
        "solverName": "potentialFlowSolver",
        "prefactor": 1.0,
        "dirichletBoundaryConditions": potential_flow_bc,
      },
    },
    "Activation": {
      "FiniteElementMethod" : {  
        "meshName": "mesh",
        "solverName": "activationSolver",
        "prefactor": 1.0,
      },
    },
    "OutputWriter" : [
      {"format": "Paraview", "outputInterval": 1, "filename": "out/output", "binary": True, "fixedFormat": False},
      #{"format": "ExFile", "filename": "out/fibre_"+str(i), "outputInterval": 1./dt_1D*output_timestep, "sphereSize": "0.02*0.02*0.02"},
      #{"format": "PythonFile", "filename": "out/fibre_"+str(i), "outputInterval": int(1./dt_1D*output_timestep), "binary":True, "onlyNodalValues":True},
    ]
  },
}
