# AMI model
A two-muscles-one-tendon model using OpenDiHu and preCICE

### How to build?
Follow opendihu's documentation for installation, then run 
```
mkorn && sr
```
### How to run?
You will need one terminal per each participant. All terminals should be at the same directory `ami/build_release`

terminal 1: left muscle
```
./muscle ../settings_muscle_left.py
```
terminal 2: tendon
```
./tendon ../settings_tendon.py
```
terminal 3: right muscle
```
./only_mechanics_muscle ../settings_muscle_right.py
```

Alternatively you can use a tendon solver that assumes a linear material: Choose `tendon_material = "SaintVenantKirchoff"` in `variables.py` and then use the executable `./linear_tendon` when running the tendon participant.

### How to run in a cluster?
Here's an example on how to lunch the three participants in the cluster (Add this code to your bash file):
```
echo "Launching left muscle"
mpirun -n 1 ./muscle_precice ../settings_muscle_left.py &> left.log &

echo "Launching tendon"
mpirun -n 1 ./tendon_precice ../settings_tendon.py &> tendon.log &

echo "Launching right muscle"
mpirun -n 1 ./muscle_mechanics_precice ../settings_muscle_right.py &> right.log

echo "Simulation completed."

```
