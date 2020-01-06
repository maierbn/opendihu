#include "interfaces/discretizable_in_time.h"

DiscretizableInTime::DiscretizableInTime()
{
}

int DiscretizableInTime::nComponentsNode()
{
  return 1;
}

void DiscretizableInTime::getComponentNames(std::vector<std::string> &componentNames)
{
  // no special component names here, this is e.g. overloaded in cellML adapter where component names are availableq
}

void DiscretizableInTime::prepareForGetOutputConnectorData()
{

}
