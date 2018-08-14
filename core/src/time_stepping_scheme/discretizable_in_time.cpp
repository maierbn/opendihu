#include "time_stepping_scheme/discretizable_in_time.h"

DiscretizableInTime::DiscretizableInTime(SolutionVectorMapping solutionVectorMapping)
 : solutionVectorMapping_(solutionVectorMapping)
{
}

int DiscretizableInTime::nComponentsNode()
{
  return 1;
}

bool DiscretizableInTime::setInitialValues(Vec& initialValues)
{
  return false;
}

SolutionVectorMapping& DiscretizableInTime::solutionVectorMapping()
{
  return solutionVectorMapping_;
}

void DiscretizableInTime::getComponentNames(std::vector<std::string> &componentNames)
{
  // no special component names here, this is e.g. overloaded in cellML adapter where component names are available
}