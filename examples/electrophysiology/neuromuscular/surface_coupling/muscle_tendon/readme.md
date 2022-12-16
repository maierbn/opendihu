# Convergence criterium
I observe that the values of displacement and traction in the interface have a wide range of numbers (from 1e-10 to 1e-2). Therefore it seems convenient to choose a relative convergence criterium instead of an absolute criterium. 

```
      <relative-convergence-measure limit="1e-2" data="Displacement" mesh="TendonMeshLeft" strict="1"/>
      <relative-convergence-measure limit="1e-2" data="Traction" mesh="MuscleMeshLeft" strict="1"/>
```

Choosing a relative criterium requires having lots of iterations for `t=0.0`.

- **no acceleration**: 
I run the simulation using `max-iterations: 250`. 
207 iterations were needed in order to converge the first time step. After that 2-5 iterations are needed per timestep, but at 10.5 we need many iterations again. In total we have 869 timesteps.  

> **Warning**
> If we look at the fibers mesh, we see that it looks like we are advancing in time during the coupling iterations. At `t=0.0` the muscle is already contracted, and the activation of the fiber fades out at `t=0.6`.

If we make the convergence criterium non-strict and we set `max-iterations: 2` and the results with the fibers look consistent with a double speed of the propagation of the action-potential. Except for the begining, the steps are generally converged.  


