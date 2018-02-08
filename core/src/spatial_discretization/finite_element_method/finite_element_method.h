#pragma once

#include "equation/type_traits.h"

#include "spatial_discretization/spatial_discretization.h"
#include "spatial_discretization/finite_element_method/02_stiffness_matrix.h"
#include "spatial_discretization/finite_element_method/04_rhs.h"
#include "spatial_discretization/finite_element_method/05_timestepping.h"
#include "basis_on_mesh/05_basis_on_mesh.h"
#include "basis_on_mesh/mixed_basis_on_mesh.h"
#include "basis_function/mixed.h"

namespace SpatialDiscretization
{
 
/** inherited class that has additional Term template parameter
 */
template<typename MeshType, typename BasisFunctionType, typename IntegratorType, typename Term, typename = Term>
class FiniteElementMethod :
  public FiniteElementMethodStiffnessMatrix<BasisOnMesh::BasisOnMesh<MeshType, BasisFunctionType>, IntegratorType, Term>
{
};

/** partial specialisation for Equation::Static::Laplace: has only stiffnessMatrix
 * use inheritage hierarchy until file 02_stiffness_matrix.h
 */
template<typename MeshType, typename BasisFunctionType, typename IntegratorType>
class FiniteElementMethod<MeshType, BasisFunctionType, IntegratorType, Equation::Static::Laplace> :
  public FiniteElementMethodStiffnessMatrix<BasisOnMesh::BasisOnMesh<MeshType, BasisFunctionType>, IntegratorType, Equation::Static::Laplace>
{
public:
  //! use constructor of base class
  using FiniteElementMethodStiffnessMatrix<BasisOnMesh::BasisOnMesh<MeshType, BasisFunctionType>, IntegratorType, Equation::Static::Laplace>
    ::FiniteElementMethodStiffnessMatrix;
 
private:
  //! initialize rhs vector to 0
  void setRightHandSide();
};

/** common class for not specialized MeshType, BasisFunctionType, for poisson equation/everything that is static and has a rhs
 * use inheritage hierarchy until file 04_rhs.h
 */
template<typename MeshType, typename BasisFunctionType, typename IntegratorType, typename Term>
class FiniteElementMethod<MeshType, BasisFunctionType, IntegratorType, Term, Equation::hasLaplaceOperatorWithRhs<Term>> :
  public FiniteElementMethodRhs<BasisOnMesh::BasisOnMesh<MeshType, BasisFunctionType>, IntegratorType, Term>
{
public:
  //! use constructor of base class
  using FiniteElementMethodRhs<BasisOnMesh::BasisOnMesh<MeshType, BasisFunctionType>, IntegratorType, Term>::FiniteElementMethodRhs;
  
};

/* class for mixed formulation for structural mechanics
 */
template<typename MeshType, typename LowOrderBasisFunctionType, typename HighOrderBasisFunctionType, typename MixedIntegratorType>
class FiniteElementMethod<MeshType, BasisFunction::Mixed<LowOrderBasisFunctionType, HighOrderBasisFunctionType>, MixedIntegratorType, Equation::Static::SolidMechanics> :
  public FiniteElementMethodRhs<BasisOnMesh::Mixed<
    BasisOnMesh::BasisOnMesh<MeshType, LowOrderBasisFunctionType>,
    BasisOnMesh::BasisOnMesh<MeshType, HighOrderBasisFunctionType>>, MixedIntegratorType, Equation::Static::SolidMechanics>
{
public:
  //! use constructor of base class
  using FiniteElementMethodRhs<BasisOnMesh::Mixed<
    BasisOnMesh::BasisOnMesh<MeshType, LowOrderBasisFunctionType>,
    BasisOnMesh::BasisOnMesh<MeshType, HighOrderBasisFunctionType>>, MixedIntegratorType, Equation::Static::SolidMechanics>::FiniteElementMethodRhs;
    
};

/** common class for not specialized MeshType, BasisFunctionType, for time stepping
 * use inheritage hierarchy until file 05_timestepping.h
 */
template<typename MeshType, typename BasisFunctionType, typename IntegratorType, typename Term>
class FiniteElementMethod<MeshType, BasisFunctionType, IntegratorType, Term, Equation::hasLaplaceOperatorWithTimeStepping<Term>> :
  public FiniteElementMethodTimeStepping<BasisOnMesh::BasisOnMesh<MeshType, BasisFunctionType>, IntegratorType, Term>
{
public:
  //! use constructor of base class
  using FiniteElementMethodTimeStepping<BasisOnMesh::BasisOnMesh<MeshType, BasisFunctionType>, IntegratorType, Term>::FiniteElementMethodTimeStepping;
  
};

}  // namespace

#include "spatial_discretization/finite_element_method/finite_element_method.tpp"
