#pragma once

#include "time_stepping_scheme/time_stepping_scheme.h"
#include "control/runnable.h"

namespace TimeSteppingScheme
{

template<typename DiscretizableInTime>
class CrankNicholson : public TimeSteppingScheme, public Runnable
{
public:
private:
};

}  // namespace