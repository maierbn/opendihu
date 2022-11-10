# Surface Coupling Examples

This folder contains examples were muscles and tendons are coupled using preCICE, a coupling library for partitioned multi-physics simulations. 
PreCICE does not provide an official adapter for OpenDiHu (yet). Benjamin Maier developed an opendihu precice adapter which is included in the opendihu source code. The existing opendihu consists of two independent adapters; one for surface coupling and another one for volume coupling.  

The examples in this folder use the built-in opendihu precice adapter for surface coupling.

## The muscle participant:

**A multi-scale problem**

The muscle participant requires of the solution of a multi-scale problem consisting of three main phenomena:
* sub-cellular processes at the sarcomeres: modelled as a 0D process and with `dt_0D = 0.5e-3` [ms]                      
* potential propagation at the fibers: modelled as a 1D process and with `dt_1D = 1e-3` [ms] 
* solid-mechanics (muscle contraction): modelled as a 3D process and with `dt_elasticity = 1e-1` [ms] 

**About muscle fibers**

The muscle fibers are created in `helper.py`.  All fibers have the same length (equal to the muscle length) and direction (given by z-axis). In the `variables.py` we can specify:
- numer of points per fiber (eg. number of sarcomeres per fiber): `n_points_whole_fiber`
- number of fibers in the x and y axis: `n_fibers_x` and `n_fibers_y`
- fiber distribution file: use input file `/MU_fibre_distribution_multidomain_67x67_100.txt`

The fibers are stimulated via cortical input. In these examples the neural estimulation is not modelled explicitely and a file with the firing times is used instead. 

- fiber firing times: use input file `"/MU_firing_times_real.txt"`
- alterative fiber firing times (zero everywhere, for no activation): use input file `"/MU_firing_times_real_no_firing.txt"`
- start time for firing: `activation_start_time`

> **Note**
> If you choose not to activate the muscle fibers the muscle will not contract unless you apply an external force. If you try this in the muscle-tendon example you will see that both the muscle and the tendon remain unmodified.

**About the subcellullar processes**

We are using the Hodgkin-Huxley-Razumova model. 

If you want to change the used cellml model you need to do the following:
- select another file for `cellml_file`
- provide the required parameters in `muscle_material_parameters`
- modify `Ǹ_states` and `Ǹ_algebraics` in `muscle_neuromuscular.cpp` according to the model

> **Warning**
> Changing the subcellular model typically requires recompiling the muscle solver.

**The structure of the muscle solver**

The opendihu solver that is used for the muscle participant is given by `muscle_neuromuscular.cpp`.

The muscle solver consists of the `MonodomainSolver`, which models the 1d voltage propagation in the muscle fibers, and the `MuscleContractionSolver`, which models the 3d solid mechanics. 

TODO: Connected slots between  `MonodomainSolver` and `MuscleContractionSolver`
The `MuscleContractionSolver` receives the following slots:
`slotNames":                    ["m1lda", "m1ldot", "m1g_in", "m1T", "m1ux", "m1uy", "m1uz"]`


The idea behind the `MonodomainSolver` is to apply strang splitting to each fiber. The components of the strang splitting are the fiber reaction term and the fiber diffusion term. The fiber reaction term corresponds to the sub-cellular processes taking place at the sarcomeres. 

> **Note**
> The overall timestep of the strang splitting is `dt_splitting_0D1D = dt_1D`
> The optimal choice of timesteps for the strang splitting is `dt_1D = 2*dt_0D`

TODO: **Connected slots between reaction term and diffusion term**

**Meshes and spatial discretization**

On the one hand we have a 1D mesh for each fiber. On the other hand we have a 3D mesh for the solid mechanics problem. The dimensions and numer of elements are defined by `muscle1_extent` ([cm]) and `n_elements_muscle1`.

The spatial discretization for the solid mechanics is done via FEM. In OpenDiHu we can choose between using linear or quadratic ansatz functions. Quadratic ansatz functions are used for the muscle participant. This is implemented by choosing `"meshName": "muscle1Mesh_quadratic"` in the muscle settings file. 

> **Note**
> Quadratic elements are created by taking 2x2x2 linear elements. Thus `n_elements_muscle1` must contain even numbers.

TODO: check that defining `"muscle1Mesh"` is actually necessary.



## The tendon participant:
Two opendihu solvers for the tendon participant can be found in this folder: `tendon_precice_dynamic.cpp` and `tendon_linear_precice_dynamic.cpp`. The main difference is that the equation describing the solid mechanics used on each of them is different. I recommend using `tendon_linear_precice_dynamic.cpp`, which solves a linear equation and is believed to be more stable. 

In any case, you must specify the `tendon_material` according to the opendihu solver you are using:

- linear tendon: `tendon_linear_precice_dynamic.cpp` and `tendon_material= "SaintVenantKirchoff"` 
- non-linear tendon: `tendon_precice_dynamic.cpp` and `tendon_material= "SaintVenantKirchoff"` 

In the file `settings_tendon.py` you can see the input parameters that are used in each case. You can play with the elasticity of the tendon by modifying these parameters. 

> **Note**
> If we solve a linear equation does it makes sense to use quadratic ansatz functions?? TODO: make sure that linear functions work well for a single tendon. 



## About the muscle-tendon example:

**How to build and run**

To build the muscle and tendon solvers:
```
cd muscle_tendon/
mkorn && sr
cd build_release
```

You will need two terminals to run the example. In the first terminal execute
```
cd muscle_tendon/build_release
./muscle_neuromuscular.cpp ../settings_muscle.py
```

and in the second one run
```
cd muscle_tendon/build_release
./tendon_linear.cpp ../settings_tendon.py
```

> **Note**
> To obtain a more detailed logging output (eg. to debug precice) it's convenient to build using debug mode. In that case run `mkorn && sd` and run the code in the `build_debug/` folder.


**Set-up**
- muscle: fixed on one end (z=0.0) and attached to the tendon in the other end (z= muscle initial length)
- tendon: fixed on one end (z= muscle initial length + tendon length) and attached to the muscle on the other end (z= muscle initial length)
- the spatial discretization of the muscle and tendon is chosen so that we have a matching mesh (we define same number of elements on the x and y direction)
- transfer of data: the tendon sends the displacement and velocities to the muscle and the muscle sends traction data to the tendon.
- `</coupling-scheme:serial-explicit>`

If you run the simulation you will see this looks quite good! 

> **Note**
> There are only some mismatches for the traction values at the edges of the coupling interface, but I believe this is due to the implementation of von Neumann boundary conditions in openDiHu and is not a worrying issue.

**Open Issues**

TODO: currently trying to get `</coupling-scheme:parallel-implicit>` to work properly

- Implicit coupling reaches the maximum number of iterations, even if it is high (eg. 100)
    - No improvement observed if acceleration schemes are used.
    - Maybe it helps using absolute convergence criterium instead of relative criterium in case this is due to an almost-zero denominator.
- As a consequence of the previous point we cannot have running simulations where the traction is sent from the tendon to the muscle.
- If we replace the free end of the tendon by a traction bc this boundary condition is not reflected in the results. However, a single tendon with different traction boundary conditions at the ends was simulated without issues.


## About the muscle-tendon-muscle example:

**How to build and run**

To build the muscle and tendon solvers:
```
cd muscle_tendon/
mkorn && sr
cd build_release
```

You will need three terminals to run the example. In the first one execute
```
cd muscle_tendon_muscle/build_release
./muscle_neuromuscular.cpp ../settings_muscle.py
```

and in the second terminal run
```
cd muscle_tendon_muscle/build_release
./tendon_linear.cpp ../settings_tendon.py
```
and in the third one run
```
cd muscle_tendon_muscle/build_release
./muscle_neuromuscular.cpp ../settings_muscle_right.py
```

**How to visualize**

The visualization can be done with *ParaView*.
Interesting group files to visualize:
- muscle1_contraction_..vtp
- muscle1_fibers_..vtp (use Points Gaussian)
- mechanics_3D_..vtp (tendon results)

**Set-up**
- Muscle left and muscle right are identical except for their location on the z axis. The two muscles are connected by a tendon.
- muscle left: fixed on one end (z=0.0) and attached to the tendon in the other end (z= muscle initial length)
- tendon: attached to muscle 2 (z= muscle initial length + tendon length) and attached to the muscle 1 on the other end (z= muscle initial length)
- muscle right: fixed on one end (z= 2* muscle initial length + tendon length) and attached to the tendon in the other end (muscle initial length + tendon length)
- the spatial discretization of the muscles and tendon is chosen so that we have a matching mesh (we define same number of elements on the x and y direction)
- transfer of data: the tendon sends displacements and velocities to the muscle and both muscles send tractions to the tendon.
- coupling scheme: 
    -   multi-coupling (this uses a combination of parallel-implicit in practice)
    -   combination of two `</coupling-scheme:serial-explicit>`


**Open Issues**
* When only the left muscle is activated:
    * the interface between the left muscle and the tendon does not move, meaning that the tendon and the right muscle remain unmodified.
* When both muscles are activated:
    * the tendon is displaced to the right. It looks like we have introduced some artificial unsimmetry to the problem. Changing the order in the xml configuration file does not change the result.
    * Using multi-coupling instead of the combination of two explicit results in the tendon being further displaced to the right. 
