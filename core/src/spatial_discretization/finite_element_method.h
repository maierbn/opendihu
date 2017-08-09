#pragma once

#include "spatial_discretization/spatial_discretization.h"
#include "control/runnable.h"
#include "control/dihu_context.h"
#include "control/problem_data.h"

namespace SpatialDiscretization
{

template<typename Mesh, typename BasisFunction, typename DiscretizationInSpace>
class FiniteElementMethod : public SpatialDiscretization, public Runnable
{
public:
  FiniteElementMethod(DihuContext &settings);
  
  // perform computation
  void run();
  static void run(ProblemData &data);
private:
  DihuContext &context_;    ///< the context object containing everything to be stored
  ProblemData data_;
};

}  // namespace

#include "spatial_discretization/finite_element_method.tpp"