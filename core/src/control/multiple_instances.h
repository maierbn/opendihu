#pragma once

#include <petscvec.h>

#include "control/runnable.h"
#include "control/dihu_context.h"
#include "data_management/solution_vector_mapping.h"

namespace Control
{

/** This class hold multiple instances of the template type, e.g. for having multiple fibres, which are each as in example electrophysiology
  */
template<class TimeSteppingScheme>
class MultipleInstances: public Runnable
{
public:
  
  //! constructor
  MultipleInstances(DihuContext context); 
 
  //! advance simulation by the given time span [startTime_, endTime_]
  void advanceTimeSpan();
 
  //! set a new time interval that will be simulated by next call to advanceTimeSpan. This also potentially changes the time step width (it preserves the number of timesteps in the new time span)
  void setTimeSpan(double startTime, double endTime);

  //! initialize time span from specificSettings_
  void initialize();
  
  //! return whether the scheme has a specified mesh type and is not independent of the mesh type
  bool knowsMeshType();
  
  //! returns the Petsc solution vector
  Vec &solution();
  
  //! return the solution vector mapping object, that contains information on if there are more internal values stored in the data_ object than may be needed for further computationo
  SolutionVectorMapping &solutionVectorMapping();

  //! run solution process
  void run();
  
protected:
  
  DihuContext context_; ///< the context object that holds the config for this class
  int nInstances_; ///< number of instances
  
  std::vector<TimeSteppingScheme> instances_;   ///< the instances of the problem
  
};

};

#include "control/multiple_instances.tpp"