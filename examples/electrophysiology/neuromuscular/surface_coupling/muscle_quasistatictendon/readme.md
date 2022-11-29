# Coupling a muscle with a quasistatic tendon

The goal of this case is to debug the use of a tendon quasistatic solver for coupling purposes. The goal is to get matching traction values at the interface, eg. get implicit coupling to converge. 

> **Warning**
> The `run()` method in the `HyperelasticitySolver` outputs after calling the non-linear solver. In order not to get the output for the tendon quasistatic duplicated you can comment the output after the solver call in `02_hyperelasticity_solver.tpp`. However I didn't push this because otherwhise we loose the output in `cubic_tendon/`.
```
  // if (withOutputWritersEnabled)
  // {
  //   this->outputWriterManager_.writeOutput(this->data_, 0, 0.0, 0);
  //   this->outputWriterManagerPressure_.writeOutput(this->pressureDataCopy_, 0, 0.0, 0);
  // }

  this->nonlinearSolve();
  postprocessSolution();

  // stop duration measurement
  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::stop(this->durationLogKey_);

  // write current output values
  // if (withOutputWritersEnabled)
  // {
  //   this->outputWriterManager_.writeOutput(this->data_, 1, endTime_);
  //   this->outputWriterManagerPressure_.writeOutput(this->pressureDataCopy_, 1, endTime_);
  // }
```

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

Assume `n_elements_muscle1 = [2, 2, 20]` and `n_elements_tendon = [2, 2, 4]`

- linear_quasistatic with `dt=0.01`:
Not looking good. The corners (xmax, 0) and (0, ymax) are off. But maybe a absolute convergence of 1e-2 could still be satisfied. 