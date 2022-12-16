# Coupling a muscle with a quasistatic tendon

The goal of this case is to debug the use of a tendon quasistatic solver for coupling purposes. The goal is to get matching traction values at the interface, eg. get implicit coupling to converge. 

## How to build and run

To build the muscle and quasistatic tendon solvers:
```
cd muscle_quasistatictendon/
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
./tendon_quasistatic.cpp ../settings_tendon.py
```

Alternatively you can choose a linear quasistatic solver for the tendon
```
cd muscle_tendon/build_release
./tendon_quasistatic.cpp ../settings_tendon.py
```

## Use of acceleration schemes
The following is considered:

```
 <max-iterations value="50" />
      <absolute-convergence-measure limit="1e-8" data="Displacement" mesh="TendonMeshLeft" strict="1"/>
      <absolute-convergence-measure limit="1e-8" data="Velocity" mesh="TendonMeshLeft" strict="1"/>
      <absolute-convergence-measure limit="1e-4" data="Traction" mesh="MuscleMeshLeft" strict="1"/>
```

> **Note**
>  t=0.7 is the first timestep were the tendon writes non-zero values of displacement and velocities. If the initial displacement is 0s everywhere, we expect the velocity to be a factor of 10 larger than the displacement
> 

- **no acceleration**: Crashes at t=0.9
- **constant acceleration**: 
| relaxation |  crashes at t |
|---|---|
| 0.1  | 0.7  | 
| 0.1 |  0.7 |  


- **aitken acceleration**: Error goes down and up.

|  initial relaxation |  crashes at t |
|---|---|
| 0.5  | 0.8  | 
| 0.1 |  0.8 |  

- **quasi-newton acceleration**:

|  initial relaxation |  preconditioner |  max-used.iterations | time-windows-reused | crashes at t | Note 
|---|---|---|---| ---|
| 0.1 | residual-sum | 8 | 20 | 0.8 | - |


1) The number of columns in the least squares system exceeded half the number of unknowns at the interface. The system will probably become bad or ill-conditioned and the quasi-Newton acceleration may not converge. Maybe the number of allowed columns ("max-used-iterations") should be limited -> since we have 9 values at the interface chance from 100 to 8. 

## About the quasistatic solver

- Add new folder with solver files `solid_mechanics/quasistatic_hyperelasticity`
    - `quasistatic_hyperelasticity_solver.h` is identical to `dynamic_hyperelasticity_solver.h`: only the solver/class name changes
    - `quasistatic_hyperelasticity_solver.tpp` is almost identical to `dynamic_hyperelasticity_solver.tpp`: in function `advanceTimeSpan(withOutputWritersEnabled)` we call function `solveQuasistaticProblem` instead of `solveDynamicProblem`
- Write function `solveQuasistaticProblem` in file `solid_mechanics/hyperelasticity/01_material_computations_wrappers.tpp`
    - add new function to file `solid_mechanics/hyperelasticity/01_material_computations_wrappers.tpp`
    - the main idea is to use the results of the previous timesteps for the initial values of displacement, velocities and pressure.
- Add precice adapter for quasistatic solver in `00_nested_solver.h` and `00_nested:solver.tpp`.

