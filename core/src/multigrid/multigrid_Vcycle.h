#pragma once

#include "operator_splitting/operator_splitting.h"

namespace Multigrid
{

template<typename FiniteElement1, typename FiniteElement2>
class multigrid_Vcycle :
  public OperatorSplitting<TimeStepping1,TimeStepping2>
{
public:
  //! constructor
  multigrid_Vcycle(DihuContext context);

  //! advance time stepping by span
  void advanceTimeSpan();

protected:
};

}  // namespace

#include "operator_splitting/multigrid_Vcycle.tpp"
