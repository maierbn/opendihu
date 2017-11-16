#pragma once

#include "time_stepping_scheme/time_stepping_scheme.h"

namespace OperatorSplitting
{

class OperatorSplitting : 
  public ::TimeSteppingScheme::TimeSteppingScheme
{
public:
  //! constructor
  OperatorSplitting(const DihuContext& context);
 
  virtual ~OperatorSplitting() {}
private:
};

}  // namespace