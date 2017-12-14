#pragma once

#include "equation/type_traits.h"

#include "spatial_discretization/spatial_discretization.h"
#include "spatial_discretization/finite_element_method/stiffness_matrix.h"
#include "spatial_discretization/finite_element_method/rhs.h"
#include "spatial_discretization/finite_element_method/timestepping.h"

namespace SpatialDiscretization
{
 
// inherited class that has additional Term template parameter
template<typename MeshType, typename BasisFunctionType, typename IntegratorType, typename Term, typename DummyForTypeTraits = Term>
class FiniteElementMethod :
  public FiniteElementMethodStiffnessMatrix<MeshType, BasisFunctionType, IntegratorType>
{
};

// partial specialisation for Equation::Static::Laplace
template<typename MeshType, typename BasisFunctionType, typename IntegratorType>
class FiniteElementMethod<MeshType, BasisFunctionType, IntegratorType, Equation::Static::Laplace> :
  public FiniteElementMethodStiffnessMatrix<MeshType, BasisFunctionType, IntegratorType>
{
public:
  using FiniteElementMethodStiffnessMatrix<MeshType, BasisFunctionType, IntegratorType>::FiniteElementMethodStiffnessMatrix;
 
private:
  void setRightHandSide();
};

// common class for not specialized MeshType, BasisFunctionType, for poisson equation
template<typename MeshType, typename BasisFunctionType, typename IntegratorType, typename Term>
class FiniteElementMethod<MeshType, BasisFunctionType, IntegratorType, Term, Equation::hasLaplaceOperatorWithRhs<Term>> :
  public FiniteElementMethodRhs<MeshType, BasisFunctionType, IntegratorType>
{
public:
  using FiniteElementMethodRhs<MeshType, BasisFunctionType, IntegratorType>::FiniteElementMethodRhs;
  
};

// common class for not specialized MeshType, BasisFunctionType, for time stepping
template<typename MeshType, typename BasisFunctionType, typename IntegratorType, typename Term>
class FiniteElementMethod<MeshType, BasisFunctionType, IntegratorType, Term, Equation::hasLaplaceOperatorWithTimeStepping<Term>> :
  public FiniteElementMethodTimeStepping<MeshType, BasisFunctionType, IntegratorType>
{
public:
  using FiniteElementMethodTimeStepping<MeshType, BasisFunctionType, IntegratorType>::FiniteElementMethodTimeStepping;
  
};


}  // namespace

#include "spatial_discretization/finite_element_method/base.tpp"
#include "spatial_discretization/finite_element_method/laplace.tpp"
#include "spatial_discretization/finite_element_method/right_hand_side.tpp"
#include "spatial_discretization/finite_element_method/timestepping.tpp"
#include "spatial_discretization/finite_element_method/stiffness_matrix_regular_fixed.tpp"
#include "spatial_discretization/finite_element_method/stiffness_matrix_deformable.tpp"