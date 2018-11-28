#pragma once

#include "operator_splitting/operator_splitting.h"

namespace OperatorSplitting
{

template<typename TimeStepping1, typename TimeStepping2>
class Strang :
  public OperatorSplitting<TimeStepping1, TimeStepping2>
{
public:
  typedef typename TimeStepping1::FunctionSpace FunctionSpace;
  typedef typename TimeStepping1::Data Data;
  
  //! constructor
  Strang(DihuContext context);

  //! advance time stepping by span
  void advanceTimeSpan();
};

}  // namespace

#include "operator_splitting/strang.tpp"
