#include "model_order_reduction/pod.h"

namespace ModelOrderReduction
{

 
template<typename DiscretizableInTimeType>
PODBase<DiscretizableInTimeType>::
PODBase(DihuContext context) : 
  DiscretizableInTime(SolutionVectorMapping(true)), problem_(context), context_(context)
{
}
 
template<typename DiscretizableInTimeType>
void PODBase<DiscretizableInTimeType>::
initialize()
{
  problem_.initialize();
}
  
//! timestepping rhs function f of equation u_t = f(u,t)
template<typename DiscretizableInTimeType>
void PODBase<DiscretizableInTimeType>::
evaluateTimesteppingRightHandSide(Vec &input, Vec &output, int timeStepNo, double currentTime)
{
  // call the method of the underlying problem
  problem_.evaluateTimesteppingRightHandSide(input, output, timeStepNo, currentTime);
}
 
template<typename DiscretizableInTimeType>
void POD<DiscretizableInTimeType, LinearPart>::
evaluateTimesteppingRightHandSide(Vec &input, Vec &output, int timeStepNo, double currentTime)
{
  PODBase<DiscretizableInTimeType>::evaluateTimesteppingRightHandSide(input, output, timeStepNo, currentTime);
};
 
//! get the number of degrees of freedom per node which is 1 by default
template<typename DiscretizableInTimeType>
int PODBase<DiscretizableInTimeType>::
nComponentsNode()
{
  return problem_.nComponentsNode();
}
  
//! set initial values and return true or don't do anything and return false
template<typename DiscretizableInTimeType>
bool PODBase<DiscretizableInTimeType>::
setInitialValues(Vec &initialValues)
{
  problem_.setInitialValues(initialValues);
}
  
//! return whether the object has a specified mesh type and is not independent of the mesh type
template<typename DiscretizableInTimeType>
bool PODBase<DiscretizableInTimeType>::
knowsMeshType()
{
  return problem_.knowsMeshType();
}
template<typename DiscretizableInTimeType>
std::shared_ptr<Mesh::Mesh> PODBase<DiscretizableInTimeType>::
mesh()
{
  return problem_.mesh();
}

};
