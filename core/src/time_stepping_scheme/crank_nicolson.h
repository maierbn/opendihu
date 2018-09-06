#pragma once

#include "time_stepping_scheme/time_stepping_scheme.h"
#include "time_stepping_scheme/time_stepping_implicit.h"

namespace TimeSteppingScheme
{

template<typename DiscretizableInTimeType>
class CrankNicholson :
public TimeSteppingImplicit<DiscretizableInTimeType>
{
public:
  
  //! constructor
  CrankNicholson(DihuContext context);
  
  //! advance simulation by the given time span [startTime_, endTime_] with given numberTimeSteps, data in solution is used, afterwards new data is in solution
  void advanceTimeSpan();
  
  void initialize();
  
protected:
  
  //! precomputes the system matrix A=(I - dt/2 * M^(-1)K) for the implicit euler scheme
  void setSystemMatrix(double timeStepWidth);

private:
  
  //! precomputes the integration matrix for the right hand side, (I + dt/2 * M^(-1)K)
  void setIntegrationMatrixRightHandSide();
  
  //! multiplies the right hand side with the integration matrix, (I + dt/2 * M^(-1)K) 
  void evaluateTimesteppingRightHandSideImplicit(Vec &input, Vec &output, int timeStepNo, double currentTime);
   
};

}  // namespace
