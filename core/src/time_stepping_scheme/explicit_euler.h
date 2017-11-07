#pragma once

#include "time_stepping_scheme/time_stepping_scheme.h"
#include "control/runnable.h"
#include "data_management/time_stepping.h"
#include "control/dihu_context.h"

namespace TimeSteppingScheme
{

template<typename DiscretizableInTime>
class ExplicitEuler : public TimeSteppingScheme, public Runnable
{
public:
 
  ///! constructor
  ExplicitEuler(DihuContext &context);
  
  ///! run the simulation
  void run();
private:
 
  void run(PyObject *specificSettings);
  
  ///! read initial values from settings and set field accordingly
  void setInitialValues(PyObject *specificSettings);
 
  DihuContext &context_;    ///< the context object containing everything to be stored
  Data::TimeStepping data_;     ///< data object that holds all PETSc vectors and matrices
  
  DiscretizableInTime discretizableInTime;    ///< the object to be discretized
  
  double endTime_;          ///< end time of simulation
  int numberTimeSteps_;     ///< number of time steps in simulation time
  double timeStepWidth_;    ///< computed time step width, endTime_/numberTimeSteps_
};

}  // namespace

#include "time_stepping_scheme/explicit_euler.tpp"