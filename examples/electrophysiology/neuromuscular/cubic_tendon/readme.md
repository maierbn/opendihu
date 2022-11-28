## Choose solver

We can choose between a dynamic or a quasistatic solver, or a linear or a non-linear solver. All in all, there are four tested solver options for this case:
- tendon_dynamic
- tendon_linear_dynamic
- tendon_quasistatic
- tendon_linear_quasistatic

Linear or non-linear is given by the material properties of the tendon and refers to the relationship between the stress S and the strain E. The material parameters must be chosen accordingly in the `variables.py`file:

- linear tendon: `tendon_material= "SaintVenantKirchoff"` 
- non-linear tendon: `tendon_material= "SaintVenantKirchoff"` 

Dynamic or static (quasistatic just means we solve a static problem for each timestep) has to do with wether we neglect inertia forces or not. For large inertia forces it is necessary to use the dynamic solver. 

> **Warning**
> The `HyperelasticitySolver` used to output before and after calling the non-linear solver. In order not to get the output duplicated I have commented the output before calling the non-linear solver.

## Bulk force and external forces

A bulk force is appled on the whole body ( eg. gravity), why a external force is a force applied to the surface.
The body force is included in the configuration by `"constantBodyForce": variables.constant_body_force` .

The external force is included in the configuration as a Neumann boundary condition. 

**Constant traction**

The most likely scenario is that we have a constant force pulling from one of the extreme of the tendon. This can be add in the configuration as follows:

```
"neumannBoundaryConditions":   variables.elasticity_neumann_bc,     
"divideNeumannBoundaryConditionValuesByTotalArea": False,         # if true we divide by the area
```

where

```
k = 0 # bc at  z=0
variables.elasticity_neumann_bc = [{"element": k*mx*my + j*mx + i, "constantVector": [0.0,0.0,-1.000], "face": "2-"} for j in range(my) for i in range(mx)]

```
Please note that tendons are rather stift, if you apply small forces (eg. < 1000) you will not see deformation with the bare eye. You see larger deformations if you use a linear solver.  

**Increasing traction**

```
"neumannBoundaryConditions":   variables.elasticity_neumann_bc,     
"divideNeumannBoundaryConditionValuesByTotalArea": False,         # if true we divide by the area
"updateNeumannBoundaryConditionsFunction": update_neumann_bc,       
"updateNeumannBoundaryConditionsFunctionCallInterval": 1          
```

where

```
external_force = 1000.0
k = 0 # bc at  z=0
variables.elasticity_neumann_bc = [{"element": k*mx*my + j*mx + i, "constantVector": [0,0,0.0], "face": "2-"} for j in range(my) for i in range(mx)]

def update_neumann_bc(t):
  factor = min(1, t/1)   # at t=1.0 we have F = external_force
  elasticity_neumann_bc = [{
		"element": k*mx*my + j*mx + i, 
		"constantVector": [0,0, -external_force*factor], 		# force pointing to bottom
		"face": "2-",
    "isInReferenceConfiguration": True
  } for j in range(my) for i in range(mx)]

  config = {
    "inputMeshIsGlobal": True,
    "divideNeumannBoundaryConditionValuesByTotalArea": False,            
    "neumannBoundaryConditions": elasticity_neumann_bc,
  }
  return config
```

> **Note**
> We use `n_elements_tendon = [6, 6, 4]` and `dt_elasticity = 0.1`. If we compare the external forces to the traction field we see that:
> - the non-linear shows more similar values than the linear
> - the quasistatic shows more similar values than the dynamic
> - The quasistatic is not matching (min is -499.2 instead -500) :disappointed:
> - TODO: TRY SMALLER TIMESTEP AND SEE IF QUASISTATIC RESULTS IMPROVED
> - TODO: COMPARE RESULTS FOR CONSTANT TRACTION


## How to set Dirichlet BC

Here we consider different possible options for the dirichlet BC assuming that we set a neumann traction boundary condition at one end of the tendon:

- **Option 1**: 

Fix one extreme of the tendon completely -> most intuitive to me
```
k = nz-1 #free side of the tendon
for j in range(ny):
  for i in range(nx):
    variables.elasticity_dirichlet_bc[k*nx*ny + j*nx + i] = [0.0, 0.0, 0.0, None, None, None]                
```

-  **Option 2**:
Fix u_z on the extreme of the tendon and fix one of the edges on the u_x and u_y direction along the tendon. 
```
k = nz-1 #free side of the tendon
for j in range(ny):
  for i in range(nx):
    variables.elasticity_dirichlet_bc[k*nx*ny + j*nx + i] = [None, None, 0.0, None, None, None] 
            
for k in range(nz):
    variables.elasticity_dirichlet_bc[k*nx*ny + j*nx][0] = 0.0   
    variables.elasticity_dirichlet_bc[k*nx*ny + j*nx][1] = 0.0                         
```

-  **Option 3**:

BM: Use *statically determined constraints* (only u_z fixed in the plane, u_x and u_y only at two edges and all three only at a corner point)

```
k = nz-1 #free side of the tendon
for j in range(ny):
  for i in range(nx):
    variables.elasticity_dirichlet_bc[k*nx*ny + j*nx + i] = [None, None, 0.0, None, None, None]       

for j in range(ny):
    variables.elasticity_dirichlet_bc[k*nx*ny + j*nx][0] = 0.0

for i in range(nx):
    variables.elasticity_dirichlet_bc[k*nx*ny  + i][1] = 0.0         
```

-  **Option 4**:

Ommit Dirichlet BC -> but then the problem is underconstraint?!


## Reproduce BM's simulation

- Choose the following parameters:

```
dt_elasticity = 0.001     
tendon_extent = [3.0, 3.0, 2.0]  //cm        
n_elements_tendon = [4, 4, 4] 
constant_body_force = None
tendon_material = "nonLinear"                 
```

> **Warning**
> Seems like `"fiberDirection": [0,0,1]` must be specified. If I change it too eg. `[0,1,0]`, results look different. Why would this make a difference?

- Initial state:

```
"initialValuesDisplacements":  [[0.0,0.0,0.0] for _ in range(nx * ny * nz)],     
"initialValuesVelocities":     [[0.0,0.0,0.0] for _ in range(nx * ny * nz)], 
```    

- Build and run:

```
mkorn && sr
./build_release/tendon_dynamic settings_tendon_bm.py
```

- Visualization in Paraview
    - u_z (representation: surface with edges)
    - T (material traction)_z (TODO:glyph?)

- Run from a user-given starting point

Define the initial values for displacement based on an input file. 

```
data = np.genfromtxt("u_1000.csv", delimiter=",")
initial_displacements = [list(x) for x in data[:,0:3]]
```

- Discussion of Results

We are solving a dynamic problem. We reach a steady state when the external forces and the internal stress field is in equilibrium. With dt = 0.001 it looks like by t=1.0 we are already in the steady state. 