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

## Discussion of results (EXPLICIT)

Assume `n_elements_muscle1 = [2, 2, 20]` and `n_elements_tendon = [2, 2, 4]`

- linear_quasistatic with `dt=0.01`:
Not looking good. The corners (xmax, 0) and (0, ymax) are off. But maybe a absolute convergence of 1e-2 could still be satisfied. 

- quasistatic with `dt=0.01`:
**Works properly** :tada:
Even the corner values match! The image below corresponds to t=10ms, when the muscle is contracted. The image shows the values provided by precice, that is, our boundary condition and the actual traction field in the tendon.

![image info](./quasistatic_tendon_at_10.png)


## Discussion of results (IMPLICIT)

We want to get implicit coupling to work so that we can reverse which participant sends which data.

Some comments:

- **Enforce one iteration**:
That would be equivalent to explicit coupling, but if we choose parallel we see that the precice values for traction are zero. This makes sense if we do a single iteration in parallel. 

- **Absolute convergence for traction 1e-4**:
It seems necessary to use the quasi-newton acceleration. In order to do the quasi-newton we have to remove the convergence criterium for the displacement, since it is always convergence. So far the best results (crashing at t=12.0) were obtained using:

```
<max-time value="20.0"/>           
<time-window-size value="0.1"/>   

<acceleration:IQN-ILS>
<!-- <data name="Displacement" mesh="TendonMeshLeft"/> -->
<data name="Traction" mesh="MuscleMeshLeft"/>
<preconditioner type="residual-sum"/>
<filter type="QR2" limit="1e-3"/>
<initial-relaxation value="0.5"/>
<max-used-iterations value="10"/>
<time-windows-reused value="20"/>
</acceleration:IQN-ILS>

<max-iterations value="50" />
<!-- <absolute-convergence-measure limit="1e-6" data="Displacement" mesh="TendonMeshLeft" strict="1"/> -->
<absolute-convergence-measure limit="1e-4" data="Traction" mesh="MuscleMeshLeft" strict="1"/>
```
