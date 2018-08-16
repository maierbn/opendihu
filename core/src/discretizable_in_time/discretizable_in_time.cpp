#include "discretizable_in_time/discretizable_in_time.h"

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