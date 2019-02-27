#pragma once

#include "operator_splitting/operator_splitting.h"

namespace OperatorSplitting
{

/**
  * A coupling between two time stepping schemes is the same as a Godunov operator splitting.
  * This class implements the functionality of both, the classes Control::Coupling and OperatorSplitting::Godunov inherit from this class.
  */
template<typename TimeStepping1, typename TimeStepping2>
class CouplingOrGodunov  :
  public OperatorSplitting<TimeStepping1,TimeStepping2>
{
public:
  //! constructor, name is the config key, i.e. "GodunovSplitting" for Godunov and "Coupling" for Coupling
  CouplingOrGodunov(DihuContext context, std::string name);

  //! advance time stepping by span
  void advanceTimeSpan();

protected:
};

}  // namespace

#include "operator_splitting/coupling_or_godunov.tpp"
