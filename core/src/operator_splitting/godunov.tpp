#include "operator_splitting/godunov.h"

namespace OperatorSplitting
{

template<typename TimeStepping1, typename TimeStepping2>
Godunov<TimeStepping1,TimeStepping2>::
Godunov(DihuContext context) :
  CouplingOrGodunov<TimeStepping1,TimeStepping2>(context, "GodunovSplitting")
{
}

};    // namespace
