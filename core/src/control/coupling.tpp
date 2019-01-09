#include "control/coupling.h"

namespace Control
{

template<typename TimeStepping1, typename TimeStepping2>
Coupling<TimeStepping1,TimeStepping2>::
Coupling(DihuContext context) :
  OperatorSplitting::CouplingOrGodunov<TimeStepping1,TimeStepping2>(context, "Coupling")
{
}

};    // namespace
