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

  //! reuse constructor of base class
  using OperatorSplitting<TimeStepping1,TimeStepping2>::OperatorSplitting;

  //! advance time stepping by span
  void advanceTimeSpan(bool withOutputWritersEnabled = true);

protected:
};

}  // namespace

#include "operator_splitting/coupling_or_godunov.tpp"
