#pragma once

#include "time_stepping_scheme/00_time_stepping_scheme.h"
#include "time_stepping_scheme/03_time_stepping_implicit.h"

namespace TimeSteppingScheme
{

template<typename DiscretizableInTimeType>
class CrankNicolson :
public TimeSteppingImplicit<DiscretizableInTimeType>
{
public:
  
  //! constructor
  CrankNicolson(DihuContext context);
  
  //! advance simulation by the given time span [startTime_, endTime_] with given numberTimeSteps, data in solution is used, afterwards new data is in solution
  void advanceTimeSpan();
  
  void initialize();

protected:
  
  //! initialize matrices and right hand side
  virtual void initializeWithTimeStepWidth_impl(double timeStepWidth);

  //! precomputes the system matrix A=(I - dt/2 * M^(-1)K) for the implicit euler scheme
  void setSystemMatrix(double timeStepWidth);

private:
  
  //! precomputes the integration matrix for the right hand side, (I + dt/2 * M^(-1)K)
  void setIntegrationMatrixRightHandSide();
  
  //! multiplies the right hand side with the integration matrix, (I + dt/2 * M^(-1)K) 
  void evaluateTimesteppingRightHandSideImplicit(Vec &input, Vec &output, int timeStepNo, double currentTime);
   
};

}  // namespace

#include "time_stepping_scheme/crank_nicolson.tpp"
