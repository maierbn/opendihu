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

## Discussion of results