This folder contains a cuboid muscle example (`MonodomainSolver`+ `MuscleContractionSolver`). Before starting the simulation, the muscle is pulled from one side using the `HyperelasticitySolver`. The other side of the muscle is fixed for the whole simulation.

To compile the code:

```
mkorn && sr
```

To run the code: 
```
cd build_release
./muscle_with_prestretch ../settings_muscle_with_prestretch.py variables_1N.py
```

There are several `variables.py` available that compute prestretch for different pulling forces. 

A requirement to run this case is to have set up `OPENDIHU_HOME` in your system. 