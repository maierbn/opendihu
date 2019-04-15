#pragma once

#include "operator_splitting/coupling_or_godunov.h"

namespace Control
{

/** Coupling is the same as Godunov splitting and is implemented in the class coupling_or_godunov.h
  */
template<typename TimeStepping1, typename TimeStepping2>
class Coupling :
  public OperatorSplitting::CouplingOrGodunov<TimeStepping1,TimeStepping2>
{
public:
  //! constructor
  Coupling(DihuContext context);
};

}  // namespace

#include "control/coupling.tpp"
