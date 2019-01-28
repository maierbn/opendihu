#include "interfaces/multipliable.h"

//! time of simulation, return 0, because there is no time stepping, this is needed in MultipleInstances to give a time stepping value to the output writer
double Multipliable::endTime()
{
  return 0;
}

//! number of time steps in simulation time, return -1 because there are no time steps
int Multipliable::numberTimeSteps()
{
  return -1;
}
