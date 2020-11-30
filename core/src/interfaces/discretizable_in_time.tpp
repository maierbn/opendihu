#include "interfaces/discretizable_in_time.h"

template<typename FunctionSpaceType, int nComponentsSolutionVariable>
void DiscretizableInTime<FunctionSpaceType,nComponentsSolutionVariable>::
getComponentNames(std::vector<std::string> &componentNames)
{
  // no special component names here, this is e.g. overloaded in cellML adapter where component names are available
}

template<typename FunctionSpaceType, int nComponentsSolutionVariable>
constexpr int DiscretizableInTime<FunctionSpaceType,nComponentsSolutionVariable>::
nComponents()
{
  return nComponentsSolutionVariable;
}

template<typename FunctionSpaceType, int nComponentsSolutionVariable>
void DiscretizableInTime<FunctionSpaceType,nComponentsSolutionVariable>::prepareForGetSlotConnectorData()
{

}
