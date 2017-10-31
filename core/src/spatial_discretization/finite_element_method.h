#pragma once

#include <Python.h>

#include "spatial_discretization/spatial_discretization.h"
#include "control/runnable.h"
#include "control/dihu_context.h"
#include "control/problem_data.h"
#include "equation/laplace.h"

namespace SpatialDiscretization
{
 
template<typename Mesh, typename BasisFunction>
class FiniteElementMethodBase : public SpatialDiscretization, public Runnable
{
public:
  FiniteElementMethodBase(DihuContext &context);
  
  // perform computation
  void run();
  void run(ProblemData &data, PyObject *settings);
  
protected:
 
  virtual void setRightHandSide(ProblemData &data, PyObject *specificSettings) = 0;
  void applyBoundaryConditions(ProblemData &data, PyObject *specificSettings);
  virtual void setStiffnessMatrix(ProblemData &data, PyObject *specificSettings) = 0;
  virtual void solve(ProblemData &data, PyObject *specificSettings);
  
  DihuContext &context_;    ///< the context object containing everything to be stored
  ProblemData data_;
};
 
// inherited class that has additional Term template parameter
template<typename Mesh, typename BasisFunction, typename Term>
class FiniteElementMethod : public FiniteElementMethodBase<Mesh, BasisFunction>
{
};

// partial specialisation for Equation::Static::Laplace
template<typename Mesh, typename BasisFunction>
class FiniteElementMethod<Mesh, BasisFunction, Equation::Static::Laplace> :
  public FiniteElementMethodBase<Mesh, BasisFunction>
{
public:
  FiniteElementMethod(DihuContext &context);
 
private:
  void setStiffnessMatrix(ProblemData &data, PyObject *specificSettings);
  void setRightHandSide(ProblemData &data, PyObject *specificSettings); 
  void applyBoundaryConditions(ProblemData &data, PyObject *specificSettings);
};

// partial specialisation for Equation::Static::Poisson
template<typename Mesh, typename BasisFunction>
class FiniteElementMethod<Mesh, BasisFunction, Equation::Static::Poisson> :
  public FiniteElementMethodBase<Mesh, BasisFunction>
{
public:
  FiniteElementMethod(DihuContext &context);
 
private:
  void setStiffnessMatrix(ProblemData &data, PyObject *specificSettings);
  void setRightHandSide(ProblemData &data, PyObject *specificSettings); 
  void applyBoundaryConditions(ProblemData &data, PyObject *specificSettings);
};



}  // namespace

#include "spatial_discretization/finite_element_method_regular_lagrange_laplace_1d.tpp"
#include "spatial_discretization/finite_element_method_regular_lagrange_laplace_2d.tpp"
#include "spatial_discretization/finite_element_method.tpp"