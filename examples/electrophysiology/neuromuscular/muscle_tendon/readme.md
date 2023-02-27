# Comparison of a single muscle with a muscle-tendon system 

## Single muscle 

This is the case setup we consider a muscle which is fixed to a wall on one end and free on the other one. The muscle is contracted due to the input from `/MU_firing_times_real.txt`. First of all, we start by studying the effect of different `dt_elasticity`. 

![image](Figure_1.png)

## Muscle-Tendon

At the free end of the muscle we add a tendon. The tendon is attached to the muscle on one end and is free on the other side. We expect the tendon to move as the muscle contracts. Since the tendon is adding to the mass that must be accelerated, it is expected to see that the muscle contracts less when the tendon is attached. 

### Explicit coupling

We use the following scheme

```
<coupling-scheme:serial-explicit>
  <participants second="TendonSolver" first="MuscleSolverLeft"/>
  <max-time value="5.0"/>           
  <time-window-size value="0.01"/>   
  <exchange data="Displacement"    mesh="TendonMeshLeft"    from="TendonSolver" to="MuscleSolverLeft"/>
  <exchange data="Velocity"    mesh="TendonMeshLeft"    from="TendonSolver" to="MuscleSolverLeft"/>
  <exchange data="Traction"        mesh="MuscleMeshLeft" from="MuscleSolverLeft" to="TendonSolver"/>  
</coupling-scheme:serial-explicit>  
```

and study the effect of the additional tendon mass:

- small tendon

```
tendon_extent = [3.0, 3.0, 4.0] # [cm, cm, cm]
tendon_offset = [0.0, 0.0, muscle1_extent[2]]
n_elements_tendon = [2, 2, 8] 
```

- large tendon

```
tendon_extent = [3.0, 3.0, 2.0] # [cm, cm, cm]
tendon_offset = [0.0, 0.0, muscle1_extent[2]]
n_elements_tendon = [2, 2, 4] 
```

![image](Figure_2.png)

### Implicit

> **Note**
> In order to use implicit coupling we have to make sure that the checkpoints are loaded correctly. In particular, it's necessary to load checkpoints for the fastmonodomainsolver, not only for the mechanics solver. An easy way to check weather the checkpoints are loaded correctly for the fibers, is to set `"outputOnlyConvergedTimeSteps": False` and check that the results in `muscle1_fibers*.vtu` do not change for every iteration. 


The next plot compares the results using different criteria for the implicit scheme. We consider cases where we do a constant number of iterations per timestep, vs cases where we set a convergence criterium based on the data. In the second case, is necessary to add an acceleration scheme. 

![image](Figure_3.png)

Note that *implicit* corresponds to the following
```
<acceleration:IQN-ILS>
  <data name="Displacement" mesh="TendonMeshLeft"/>
  <data name="Velocity" mesh="TendonMeshLeft"/>
  <data name="Traction" mesh="MuscleMeshLeft"/>
  <preconditioner type="residual-sum"/>
  <filter type="QR2" limit="1e-3"/>
  <initial-relaxation value="0.6"/>
  <max-used-iterations value="15"/>
  <time-windows-reused value="15"/>
</acceleration:IQN-ILS>

<max-iterations value="100"/>
<!-- <min-iteration-convergence-measure min-iterations="4" data="Displacement" mesh="TendonMeshLeft" strict="0" suffices="0"/> -->


<relative-convergence-measure limit="1e-5" data="Displacement" mesh="TendonMeshLeft" strict="1"/>
<relative-convergence-measure limit="1e-5" data="Velocity" mesh="TendonMeshLeft" strict="1"/>
<absolute-convergence-measure limit="1e-2" data="Traction" mesh="MuscleMeshLeft" strict="1"/>
```

If we look at the total number of iterations, *implicit* required 1393 iterations, which is notably less than for *implicit (4 iterations)*, where 4x499=1996 and yet the convergence measurements where not satisfied.

We did a parameter study for *implicit*:

|  `initial-relaxation value` |  # total iterations |  
|---|---|
| 0.5  | 1317  |   
| 0.6  |  1393 |   
| 0.7  |  1444 |   
| 0.8  |  1274 |   
