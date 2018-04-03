#pragma once

#include "spatial_discretization/finite_element_method/00_base.h"
#include "spatial_discretization/finite_element_method/solid_mechanics/solid_mechanics_utility.h"
#include "equation/mooney_rivlin_incompressible.h"
#include "equation/type_traits.h"

namespace SpatialDiscretization
{

/** specialisation for incompressible solid mechanics, mixed formulation with static condensation
 */
/** Helper template
 */
template<typename HighOrderBasisOnMeshType, int completePolynomialOrder>
using MixedBasisOnMeshTemplate = 
BasisOnMesh::Mixed<
  BasisOnMesh::BasisOnMesh<
    typename HighOrderBasisOnMeshType::Mesh,
    BasisFunction::CompletePolynomialOfDimensionAndOrder<HighOrderBasisOnMeshType::Mesh::dim(),completePolynomialOrder>
  >,
  HighOrderBasisOnMeshType
>;

/** specialisation for incompressible solid mechanics, mixed formulation with static condensation
 */
template<typename HighOrderBasisOnMeshType, int completePolynomialOrder, typename MixedQuadratureType, typename Term>
class FiniteElementMethodStiffnessMatrix<
  MixedBasisOnMeshTemplate<HighOrderBasisOnMeshType, completePolynomialOrder>,
  MixedQuadratureType,
  Term,
  Mesh::isDeformable<typename HighOrderBasisOnMeshType::Mesh>,
  Equation::isIncompressible<Term>
> :
  public FiniteElementMethodBase<MixedBasisOnMeshTemplate<HighOrderBasisOnMeshType, completePolynomialOrder>, MixedQuadratureType>,
  public SolidMechanicsUtility<MixedBasisOnMeshTemplate<HighOrderBasisOnMeshType, completePolynomialOrder>, MixedQuadratureType, Term>
{
public:
  typedef MixedBasisOnMeshTemplate<HighOrderBasisOnMeshType, completePolynomialOrder> MixedBasisOnMesh;
 
  // use constructor of base class
  using FiniteElementMethodBase<MixedBasisOnMesh, MixedQuadratureType>::FiniteElementMethodBase;
  
  //! assemble the stiffness matrix
  void setStiffnessMatrix();
  
  //! no set right hand side functionality
  void setRightHandSide(){}
};

/** specialisation for incompressible solid mechanics, mixed formulation without static condensation
 */
template<typename LowOrderBasisOnMeshType, typename HighOrderBasisOnMeshType, typename MixedQuadratureType, typename Term>
class FiniteElementMethodStiffnessMatrix<
  BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>,
  MixedQuadratureType, 
  Term,
  Mesh::isDeformable<typename HighOrderBasisOnMeshType::Mesh>,
  Equation::isIncompressible<Term>
> :
  public FiniteElementMethodBase<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>, MixedQuadratureType>,
  public SolidMechanicsUtility<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>, MixedQuadratureType, Term>
{
public:
  typedef BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType> MixedBasisOnMesh;
 
  // use constructor of base class
  using FiniteElementMethodBase<MixedBasisOnMesh, MixedQuadratureType>::FiniteElementMethodBase;
  
  //! assemble the stiffness matrix
  void setStiffnessMatrix();
  
  //! no set right hand side functionality
  void setRightHandSide(){}
};

/** specialisation for incompressible solid mechanics, not mixed formulation, i.e. penalty formulation,
 * currently not implemented
 *//*
template<typename BasisOnMeshType, typename QuadratureType, typename Term>
class FiniteElementMethodStiffnessMatrix<
  BasisOnMeshType,
  QuadratureType,
  Term,
  Mesh::isDeformable<typename BasisOnMeshType::Mesh>,
  Equation::isIncompressible<Term>
> :
  public FiniteElementMethodBase<BasisOnMeshType, QuadratureType>,
  public SolidMechanicsUtility<BasisOnMeshType, QuadratureType, Term>
{
public:
  // use constructor of base class
  using FiniteElementMethodBase<BasisOnMeshType, QuadratureType>::FiniteElementMethodBase;
  
  //! assemble the stiffness matrix
  void setStiffnessMatrix();
  
  //! no set right hand side functionality
  void setRightHandSide(){}
};*/

};  // namespace

#include "spatial_discretization/finite_element_method/solid_mechanics/02_stiffness_matrix_incompressible_mixed.tpp"
#include "spatial_discretization/finite_element_method/solid_mechanics/02_stiffness_matrix_incompressible_mixed_condensation.tpp"
//#include "spatial_discretization/finite_element_method/solid_mechanics/02_stiffness_matrix_incompressible_penalty.tpp"
