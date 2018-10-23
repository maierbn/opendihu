#include "interfaces/splitable.h"

Splitable::Splitable()
{
  solutionVectorMapping_ = std::make_shared<SolutionVectorMapping>();
}

std::shared_ptr<SolutionVectorMapping> Splitable::solutionVectorMapping()
{
  return solutionVectorMapping_;
}
