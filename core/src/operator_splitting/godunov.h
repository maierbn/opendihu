#pragma once

#include "operator_splitting/operator_splitting.h"

namespace OperatorSplitting
{

template<typename TimeStepping1, typename TimeStepping2>
class Godunov :
  public OperatorSplitting<TimeStepping1,TimeStepping2>
{
public:
  //! constructor
  Godunov(DihuContext context);

  //! advance time stepping by span
  void advanceTimeSpan();

protected:
};

}  // namespace

#include "operator_splitting/godunov.tpp"
