# Convergence criterium
I observe that the values of displacement and traction in the interface have a wide range of numbers (from 1e-10 to 1e-2). Therefore it seems convenient to choose a relative convergence criterium instead of an absolute criterium. 

```
      <relative-convergence-measure limit="1e-2" data="Displacement" mesh="TendonMeshLeft" strict="1"/>
      <relative-convergence-measure limit="1e-2" data="Traction" mesh="MuscleMeshLeft" strict="1"/>
```

Choosing a relative criterium makes the convergence at `t=0.0` difficult, lots of iterations are needed. 

- **no acceleration**: 
I run the simulation using `max-iterations: 300`. 207 iterations were needed in order to converge the first time step. After that 2-5 iterations are needed per timestep, but at 10.5 we need many iterations again. 

- **constant acceleration**: 
| relaxation |  crashes at t |
|---|---|
| 0.1  | 0.6  |


- **aitken acceleration**: Error goes down and up.

|  initial relaxation |  crashes at t |
|---|---|
| 0.5  | 0.8  | 
 

- **quasi-newton acceleration**:

|  initial relaxation |  preconditioner |  max-used-iterations | time-windows-reused | crashes at t | Note 
|---|---|---|---| ---|
| 0.1 | residual-sum | 8 | 100 |  | - |


1) The number of columns in the least squares system exceeded half the number of unknowns at the interface. The system will probably become bad or ill-conditioned and the quasi-Newton acceleration may not converge. Maybe the number of allowed columns ("max-used-iterations") should be limited -> since we have 9 values at the interface chance from 100 to 8. 
