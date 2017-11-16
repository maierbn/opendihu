#pragma once

#include "control/dihu_context.h"
#include "data_management/solution_vector_mapping.h"
#include "data_management/time_stepping.h"
#include "time_stepping_scheme/time_stepping_scheme.h"
#include "data_management/data.h"

namespace TimeSteppingScheme
{

template<typename DiscretizableInTime>
class TimeSteppingSchemeOde : public TimeSteppingScheme
{
public:
  TimeSteppingSchemeOde(const DihuContext &context); 
 
  //! destructor  
  virtual ~TimeSteppingSchemeOde() {}

  //! run simulation
  virtual void run();
  
  //! get the solutionVectorMapping object that stores information about which values of the solution should be used for further computation and how they can be retrieved
  SolutionVectorMapping &solutionVectorMapping();
  
  //! return the Petsc solution vector
  Vec &solution();
  
  //! return the data object
  Data::Data &data();
  
  //! initialize discretizableInTime
  void initialize();
  
protected:
  
  //! read initial values from settings and set field accordingly
  void setInitialValues();

  Data::TimeStepping data_;     ///< data object that holds all PETSc vectors and matrices
  
  int timeStepOutputFrequency_;    ///< time step number and time is output every timeStepOutputFrequency_ time steps
  DiscretizableInTime discretizableInTime;    ///< the object to be discretized
};
}  // namespace

#include "time_stepping_scheme/time_stepping_scheme_ode.tpp"
