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
  public FiniteElementMethodBase<MixedBasisOnMeshTemplate<HighOrderBasisOnMeshType, completePolynomialOrder>, MixedQuadratureType, Term>,
  public SolidMechanicsUtility<HighOrderBasisOnMeshType, Term>
{
public:
  typedef MixedBasisOnMeshTemplate<HighOrderBasisOnMeshType, completePolynomialOrder> MixedBasisOnMesh;
 
  // use constructor of base class
  using FiniteElementMethodBase<MixedBasisOnMesh, MixedQuadratureType, Term>::FiniteElementMethodBase;
  
  //! assemble the stiffness matrix
  void setStiffnessMatrix();

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
  public FiniteElementMethodBase<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>, MixedQuadratureType, Term>,
  public SolidMechanicsUtility<HighOrderBasisOnMeshType, Term>
{
public:
  typedef BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType> MixedBasisOnMesh;
 
  // use constructor of base class
  using FiniteElementMethodBase<MixedBasisOnMesh, MixedQuadratureType, Term>::FiniteElementMethodBase;
  
  //! assemble the stiffness matrix
  void setStiffnessMatrix();
};

/** specialisation for incompressible solid mechanics, not mixed formulation, i.e. penalty formulation,
 * currently not implemented
 */
template<typename BasisOnMeshType, typename QuadratureType, typename Term>
class FiniteElementMethodStiffnessMatrix<
  BasisOnMeshType,
  QuadratureType,
  Term,
  Mesh::isDeformable<typename BasisOnMeshType::Mesh>,
  Equation::isIncompressible<Term>
> :
  public FiniteElementMethodBase<BasisOnMeshType, QuadratureType, Term>,
  public SolidMechanicsUtility<BasisOnMeshType, Term>
{
public:
  // use constructor of base class
  using FiniteElementMethodBase<BasisOnMeshType, QuadratureType, Term>::FiniteElementMethodBase;
  
  //! assemble the stiffness matrix
  void setStiffnessMatrix();
  
  //! set current value of displacements, called from a PETSc SNES callback
  void setDisplacements(Vec &x);
  
  //! return the tangent stiffness matrix, called from a PETSc SNES callback
  Vec &tangentStiffnessMatrix();
  
  //! compute the internal virtual work term, dW_int
  void computeInternalVirtualWork(Vec &result);
  
  //! set entries in f to the entry in rhs for which Dirichlet BC are set
  void applyDirichletBoundaryConditionsInNonlinearFunction(Vec &f);
  
protected:
  //! solve nonlinear system
  virtual void solve();
  
  //! initialize everything, set rhs and compute tangent stiffness matrix
  virtual void initialize();
  
  //! initialize Dirichlet boundary conditions
  void initializeBoundaryConditions();
  
  std::vector<dof_no_t> dirichletIndices_;  ///< the indices of unknowns (not dofs) for which the displacement is fixed
  std::vector<double> dirichletValues_;     ///< the to dirichletIndices corresponding fixed values for the displacement
  std::vector<double> rhsValues_;           ///< the entries in rhs for which for u Dirichlet values are specified
};

};  // namespace

#include "spatial_discretization/finite_element_method/solid_mechanics/02_stiffness_matrix_incompressible_mixed.tpp"
#include "spatial_discretization/finite_element_method/solid_mechanics/02_stiffness_matrix_incompressible_mixed_condensation.tpp"
#include "spatial_discretization/finite_element_method/solid_mechanics/02_stiffness_matrix_incompressible_penalty.tpp"
#include "spatial_discretization/finite_element_method/solid_mechanics/nonlinear_solve.tpp"
