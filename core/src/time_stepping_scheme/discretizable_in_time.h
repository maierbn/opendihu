#pragma once

#include <petscsys.h>
#include <petscksp.h>

#include "data_management/solution_vector_mapping.h"

class DiscretizableInTime
{
public:
  //! constructor
  DiscretizableInTime(SolutionVectorMapping solutionVectorMapping);
 
  //! initialize timestepping
  virtual void initialize() = 0;
  
  //! timestepping rhs function f of equation u_t = f(u,t)
  virtual void evaluateTimesteppingRightHandSide(Vec &input, Vec &output, int timeStepNo, double currentTime) = 0;
 
  //! get the number of degrees of freedom per node which is 1 by default
  virtual int numberDegreesOfFreedomPerNode();
  
  //! set initial values and return true or don't do anything and return false
  virtual bool setInitialValues(Vec &initialValues);
  
  //! return the solution vector mapping object, that contains information on if there are more internal values stored in the data_ object than may be needed for further computationo
  SolutionVectorMapping &solutionVectorMapping();
  
protected:
  SolutionVectorMapping solutionVectorMapping_;   ///< the solution vector mapping object that contains information if for further computation only a subset of the stored entries in the data_.solution vector will be needed
};