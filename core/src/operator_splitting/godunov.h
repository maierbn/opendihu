#pragma once

#include "operator_splitting/operator_splitting.h"

namespace OperatorSplitting
{

/** Godunov splitting is implemented in the class coupling_or_godunov.h
  */
template<typename TimeStepping1, typename TimeStepping2>
class Godunov :
  public CouplingOrGodunov<TimeStepping1,TimeStepping2>
{
public:
  //! constructor
  Godunov(DihuContext context);
};

}  // namespace

#include "operator_splitting/godunov.tpp"
