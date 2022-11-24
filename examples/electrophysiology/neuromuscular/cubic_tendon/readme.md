## Choose solver

There are two tested solver options for this case:
- tendon_dynamic
- tendon_linear_dynamic

The tendon_linear_dynamic is less costly. For "small" deformations, the linear solver should be sufficient. 

> **Warning**
> **TODO:** Find a static solver that can run over time!

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