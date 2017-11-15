#pragma once

#include <petscsys.h>
#include <petscksp.h>

class DiscretizableInTime
{
public:
  //! initialize timestepping
  virtual void initialize() = 0;
  
  //! timestepping rhs function f of equation u_t = f(u,t)
  virtual void evaluateTimesteppingRightHandSide(Vec &input, Vec &output, int timeStepNo, double currentTime) = 0;
 
  //! get the number of degrees of freedom per node which is 1 by default
  virtual int numberDegreesOfFreedomPerNode();
  
  //! set initial values and return true or don't do anything and return false
  virtual bool setInitialValues(Vec &initialValues);
  
private:
};