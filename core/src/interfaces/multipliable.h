#pragma once

#include <Python.h>  // has to be the first included header
#include <iostream>

/**
 *  Class that can be part of a MultipleInstances
 */
class Multipliable
{
public:
  //! time of simulation, return 0, because there is no time stepping, this is needed in MultipleInstances to give a time stepping value to the output writer
  virtual double endTime();

  //! number of time steps in simulation time, return -1 because there are no time steps
  virtual int numberTimeSteps();

protected:
};
