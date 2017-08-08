#pragma once

#include "time_stepping_scheme/time_stepping_scheme.h"
#include "control/runnable.h"

namespace TimeSteppingScheme
{

template<typename DiscretizableInTime>
class ExplicitEuler : public TimeSteppingScheme, Runnable
{
public:
private:
};

}  // namespace