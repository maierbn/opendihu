#include "control/coupling/coupling.h"

namespace Control
{

template<typename TimeStepping1, typename TimeStepping2>
Coupling<TimeStepping1,TimeStepping2>::
Coupling(DihuContext context) :
  OperatorSplitting::CouplingOrGodunov<TimeStepping1,TimeStepping2>(context, "Coupling")
{
}

template<typename TimeStepping1, typename TimeStepping2>
Coupling<TimeStepping1,TimeStepping2>::
Coupling(DihuContext context, TimeStepping1 &&timeStepping1, TimeStepping2 &&timeStepping2) :
  OperatorSplitting::CouplingOrGodunov<TimeStepping1,TimeStepping2>(context, "Coupling", std::move(timeStepping1), std::move(timeStepping2))
{
}

}  // namespace
