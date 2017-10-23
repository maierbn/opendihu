#pragma once

#include <Python.h>

#include "spatial_discretization/spatial_discretization.h"
#include "control/runnable.h"
#include "control/dihu_context.h"
#include "control/problem_data.h"

namespace SpatialDiscretization
{

template<typename Mesh, typename BasisFunction, typename Term>
class FiniteElementMethod : public SpatialDiscretization, public Runnable
{
public:
  FiniteElementMethod(DihuContext &context);
  
  // perform computation
  void run();
  void run(ProblemData &data, PyObject *settings);
private:
 
  void setRightHandSide(ProblemData &data, PyObject *specificSettings); 
  void applyBoundaryConditions(ProblemData &data, PyObject *specificSettings);
  void setStiffnessMatrix(ProblemData &data, PyObject *specificSettings);
  void solve(ProblemData &data, PyObject *specificSettings);
  
  DihuContext &context_;    ///< the context object containing everything to be stored
  ProblemData data_;
};

}  // namespace

#include "spatial_discretization/finite_element_method_regular_lagrange_laplace_1d.tpp"
#include "spatial_discretization/finite_element_method_regular_lagrange_laplace_2d.tpp"
#include "spatial_discretization/finite_element_method.tpp"