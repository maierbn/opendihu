#pragma once

#include "control/load_balancing/load_balancing_base.h"
#include "operator_splitting/strang.h"
#include "time_stepping_scheme/heun_adaptiv.h"

namespace Control
{

/** This is the load balancing class for any time stepping type, it does not perform load balancing.
  */
template<typename TimeStepping>
class LoadBalancing:
  public LoadBalancingBase<TimeStepping>
{
public:

  //! use constructor of base class
  using LoadBalancingBase<TimeStepping>::LoadBalancingBase;
};

/** This is the partial specialization for strang splitting, for which load balancing is implemented.
  */
template<typename CellMLAdapter, typename DiffusionTimeStepping>
class LoadBalancing<
  OperatorSplitting::Strang<
    TimeSteppingScheme::HeunAdaptiv<CellMLAdapter>,
    DiffusionTimeStepping
  >> :
public LoadBalancingBase<OperatorSplitting::Strang<TimeSteppingScheme::HeunAdaptiv<CellMLAdapter>,DiffusionTimeStepping>>
{
public:

  //! constructor
  LoadBalancing(DihuContext context);

protected:

  //! check if the degrees of freedom should be redistributed among the ranks
  virtual void rebalance();

private:

  // Counter variable to keep track of rebalancing
  double rebalanceCounter_;

  // Frequency of the rebalancing
  int rebalanceFrequency_;
};

}  // namespace

#include "control/load_balancing/load_balancing.tpp"
