#pragma once

#include "equation/type_traits.h"

#include "spatial_discretization/spatial_discretization.h"
#include "spatial_discretization/finite_element_method/01_matrix.h"
//#include "spatial_discretization/finite_element_method/solid_mechanics/02_stiffness_matrix_compressible.h"
//#include "spatial_discretization/finite_element_method/solid_mechanics/02_stiffness_matrix_incompressible.h"
#include "spatial_discretization/finite_element_method/04_rhs.h"
#include "spatial_discretization/finite_element_method/05_time_stepping.h"
#include "function_space/function_space.h"
#include "function_space/mixed_function_space.h"
#include "basis_function/mixed.h"


#include "control/dihu_context.h"

namespace SpatialDiscretization
{

/** inherited class that has additional Term template parameter
 */
template<typename MeshType, typename BasisFunctionType, typename QuadratureType, typename Term, typename = Term, typename = BasisFunctionType>
class FiniteElementMethod :
  public FiniteElementMethodMatrix<FunctionSpace::FunctionSpace<MeshType, BasisFunctionType>, QuadratureType, Term>
{
public:
};

/** partial specialisation for Laplace: has only stiffnessMatrix
 * use inheritage hierarchy until file 03_boundary_conditions.h
 */
template<typename MeshType, typename BasisFunctionType, typename QuadratureType, typename Term>
class FiniteElementMethod<MeshType, BasisFunctionType, QuadratureType, Term, Equation::hasNoRhs<Term>, BasisFunction::isNotMixed<BasisFunctionType>> :
  public BoundaryConditions<FunctionSpace::FunctionSpace<MeshType, BasisFunctionType>, QuadratureType, Term>
{
public:
  //! use constructor of base class
  using BoundaryConditions<FunctionSpace::FunctionSpace<MeshType, BasisFunctionType>, QuadratureType, Term>::BoundaryConditions;

protected:
  //! initialize rhs vector to 0
  void setRightHandSide();
};

/* class for mixed formulation for structural mechanics
 */
template<typename MeshType, typename LowOrderBasisFunctionType, typename HighOrderBasisFunctionType, typename MixedQuadratureType, typename Term>
class FiniteElementMethod<MeshType, BasisFunction::Mixed<LowOrderBasisFunctionType, HighOrderBasisFunctionType>, MixedQuadratureType, Term> :
  public FiniteElementMethodRhs<FunctionSpace::Mixed<
    FunctionSpace::FunctionSpace<MeshType, LowOrderBasisFunctionType>,
    FunctionSpace::FunctionSpace<MeshType, HighOrderBasisFunctionType>>, MixedQuadratureType, Term>
{
public:
  //! use constructor of base class
  using FiniteElementMethodRhs<FunctionSpace::Mixed<
    FunctionSpace::FunctionSpace<MeshType, LowOrderBasisFunctionType>,
    FunctionSpace::FunctionSpace<MeshType, HighOrderBasisFunctionType>>, MixedQuadratureType, Term>::FiniteElementMethodRhs;

protected:
  void setRightHandSide(){}
};

/** common class for not specialized MeshType, BasisFunctionType, for poisson equation/solid mechanics/everything that is static and has a rhs
 * use inheritage hierarchy until file 04_rhs.h
 */
template<typename MeshType, typename BasisFunctionType, typename QuadratureType, typename Term>
class FiniteElementMethod<MeshType, BasisFunctionType, QuadratureType, Term, Equation::hasRhsNoTimestepping<Term>, BasisFunction::isNotMixed<BasisFunctionType>> :
  public FiniteElementMethodRhs<FunctionSpace::FunctionSpace<MeshType, BasisFunctionType>, QuadratureType, Term>
{
public:
  //! use constructor of base class
  using FiniteElementMethodRhs<FunctionSpace::FunctionSpace<MeshType, BasisFunctionType>, QuadratureType, Term>::FiniteElementMethodRhs;

};

/** common class for not specialized MeshType, BasisFunctionType, for time stepping
 * use inheritage hierarchy until file 05_timestepping.h
 */
template<typename MeshType, typename BasisFunctionType, typename QuadratureType, typename Term>
class FiniteElementMethod<MeshType, BasisFunctionType, QuadratureType, Term, Equation::usesTimeStepping<Term>> :
  public FiniteElementMethodTimeStepping<FunctionSpace::FunctionSpace<MeshType, BasisFunctionType>, QuadratureType, Term>
{
public:
  //! use constructor of base class
  using FiniteElementMethodTimeStepping<FunctionSpace::FunctionSpace<MeshType, BasisFunctionType>, QuadratureType, Term>::FiniteElementMethodTimeStepping;

};

}  // namespace

#include "spatial_discretization/finite_element_method/finite_element_method.tpp"
