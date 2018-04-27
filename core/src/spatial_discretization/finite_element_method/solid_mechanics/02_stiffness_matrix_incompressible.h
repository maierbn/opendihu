#pragma once

#include <vector>

#include "spatial_discretization/finite_element_method/00_base.h"
#include "spatial_discretization/finite_element_method/solid_mechanics/solid_mechanics_utility.h"
#include "spatial_discretization/finite_element_method/02_stiffness_matrix.h"
#include "equation/mooney_rivlin_incompressible.h"
#include "equation/type_traits.h"
#include "mesh/face_t.h"

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
  std::enable_if_t<LowOrderBasisOnMeshType::BasisFunction::isNodalBased, typename HighOrderBasisOnMeshType::Mesh>,
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
 */
template<typename BasisOnMeshType, typename QuadratureType, typename Term>
class FiniteElementMethodStiffnessMatrix<
  BasisOnMeshType,
  QuadratureType,
  Term,
  typename BasisOnMeshType::Mesh,
  Equation::isIncompressible<Term>,
  BasisFunction::isNotMixed<typename BasisOnMeshType::BasisFunction>
> :
  public FiniteElementMethodBase<BasisOnMeshType, QuadratureType, Term>,
  public SolidMechanicsUtility<BasisOnMeshType, Term>
{
public:
  // use constructor of base class
  using FiniteElementMethodBase<BasisOnMeshType, QuadratureType, Term>::FiniteElementMethodBase;
  
  //! assemble the stiffness matrix
  void setStiffnessMatrix(Mat stiffnessMatrix);
  void setStiffnessMatrix(){setStiffnessMatrix(PETSC_NULL);}
  
  //! return the tangent stiffness matrix, called from a PETSc SNES callback
  Mat &tangentStiffnessMatrix();
  
  //! return the virtual external energy which is constant
  Vec &rightHandSide();
  
  //! return the displacements values
  Vec &displacements();
  
  //! set the internal displacement variable as copy of the given values
  void setDisplacements(Vec &solverDisplacementVariable);
  
  //! compute the value of dW_int - dW_ext
  void computeInternalMinusExternalVirtualWork(Vec &result);
  
  //! set entries in displacements to Dirichlet BC values
  void applyDirichletBoundaryConditionsInDisplacements();
  
  //! set entries in f to the entry in rhs for which Dirichlet BC are set
  void applyDirichletBoundaryConditionsInNonlinearFunction(Vec &f);
  
  //! set rows and columns in stiffness matrix to 0 for which boundary conditions are specified
  void applyDirichletBoundaryConditionsInStiffnessMatrix(Mat &matrix);
  
  //! copy all values that are not constrained by dirichlet BC nor are z-displacements for 2D problems from the input to the output vector
  void reduceVector(Vec &input, Vec &output);
  
  //! reverse operation to reduceVector
  void expandVector(Vec &input, Vec &output);
  
  //! compute and return the appropriate analytical stiffness matrix
  void computeAnalyticalStiffnessMatrix(Mat &solverStiffnessMatrix);
  
  //! initialize everything, set rhs and compute tangent stiffness matrix
  virtual void initialize();
  
  //! print boundary conditions
  void printBoundaryConditions();  
  
protected:
  //! solve nonlinear system
  virtual void solve() override;
  
  //! debugging method
  void debug();
  
  //! initialize Dirichlet boundary conditions
  void initializeBoundaryConditions();
  
  //! read material parameters from config and set the values for static expressions within SEMT
  void initializeMaterialParameters();
  
  //! compute the internal virtual work term, dW_int
  void computeInternalVirtualWork(Vec &result);
  
  //! get the extern virtual work term, dW_ext. It is computed if it depends on displacements, otherwise the stored value is returned.
  void getExternalVirtualWork(Vec &result);
  
  //! compute the extern virtual work term, dW_ext.  
  void computeExternalVirtualWork(Vec &result);
  
  //! extract the submatrix that only contains entries for dofs that are not constraint by Dirichlet BCs
  void reduceMatrix(Mat &input, Mat &output);
  
  std::vector<dof_no_t> dirichletIndices_;  ///< the indices of unknowns (not dofs) for which the displacement is fixed
  std::vector<double> dirichletValues_;     ///< the to dirichletIndices corresponding fixed values for the displacement
  std::vector<double> zeros_;           ///< a vector of 0s, number of dirichlet values
  
  //TODO split into boundary conditions class
  struct TractionBoundaryCondition
  {
    element_no_t elementGlobalNo;
    
    Mesh::face_t face;
    std::vector<std::pair<dof_no_t, VecD<BasisOnMeshType::dim()>>> dofVectors;  //<element-local dof no, value>
    
    // parse values from python config, e.g. {"element": 1, "face": "0+", "dofVectors:", {0: [tmax,0,0], 1: [tmax,0,0], 2: [tmax,0,0], 3: [tmax,0,0]}}
    TractionBoundaryCondition(PyObject *specificSettings, std::shared_ptr<BasisOnMeshType> mesh);
  };
  
  std::vector<TractionBoundaryCondition> tractionReferenceConfiguration_;   //< tractions for elements
  std::vector<TractionBoundaryCondition> tractionCurrentConfiguration_;
  
  std::vector<std::pair<element_no_t, VecD<BasisOnMeshType::dim()>>> bodyForceReferenceConfiguration_;  //< <element global no, vector>
  std::vector<std::pair<element_no_t, VecD<BasisOnMeshType::dim()>>> bodyForceCurrentConfiguration_;    //< <element global no, vector>
  
  bool tangentStiffnessMatrixInitialized_ = false;   ///< if the tangentStiffnessMatrix has been set 
};

};  // namespace

#include "spatial_discretization/finite_element_method/solid_mechanics/02_stiffness_matrix_incompressible_mixed.tpp"
#include "spatial_discretization/finite_element_method/solid_mechanics/02_stiffness_matrix_incompressible_mixed_condensation.tpp"
#include "spatial_discretization/finite_element_method/solid_mechanics/02_stiffness_matrix_incompressible_penalty.tpp"
#include "spatial_discretization/finite_element_method/solid_mechanics/nonlinear_solve.tpp"
