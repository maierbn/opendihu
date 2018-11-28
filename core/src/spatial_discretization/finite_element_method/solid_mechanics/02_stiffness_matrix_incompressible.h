#pragma once

#include <vector>

#include "spatial_discretization/finite_element_method/00_base.h"
#include "spatial_discretization/finite_element_method/solid_mechanics/solid_mechanics_boundary_conditions.h"
#include "spatial_discretization/finite_element_method/solid_mechanics/solid_mechanics_common.h"
#include "spatial_discretization/finite_element_method/01_matrix.h"
#include "equation/mooney_rivlin_incompressible.h"
#include "equation/type_traits.h"
#include "mesh/face_t.h"

namespace SpatialDiscretization
{

/** specialisation for incompressible solid mechanics, mixed formulation with static condensation
 */
/** Helper template
 */
template<typename HighOrderFunctionSpaceType, int completePolynomialOrder>
using MixedFunctionSpaceTemplate =
FunctionSpace::Mixed<
  FunctionSpace::FunctionSpace<
    typename HighOrderFunctionSpaceType::Mesh,
    BasisFunction::CompletePolynomialOfDimensionAndOrder<HighOrderFunctionSpaceType::Mesh::dim(),completePolynomialOrder>
  >,
  HighOrderFunctionSpaceType
>;

/** u-p Mixed formulation with static condensation of the pressure.
 * Specialisation for incompressible solid mechanics, mixed formulation with static condensation
 */
template<typename HighOrderFunctionSpaceType, int completePolynomialOrder, typename MixedQuadratureType, typename Term>
class FiniteElementMethodMatrix<
  MixedFunctionSpaceTemplate<HighOrderFunctionSpaceType, completePolynomialOrder>,
  MixedQuadratureType,
  Term,
  Mesh::isDeformable<typename HighOrderFunctionSpaceType::Mesh>,
  Equation::isIncompressible<Term>
> :
  public FiniteElementMethodBase<MixedFunctionSpaceTemplate<HighOrderFunctionSpaceType, completePolynomialOrder>, MixedQuadratureType, Term>,
  public SolidMechanicsUtility<HighOrderFunctionSpaceType, Term>
{
public:
  typedef MixedFunctionSpaceTemplate<HighOrderFunctionSpaceType, completePolynomialOrder> MixedFunctionSpace;

  // use constructor of base class
  using FiniteElementMethodBase<MixedFunctionSpace, MixedQuadratureType, Term>::FiniteElementMethodBase;

  //! assemble the stiffness matrix
  void setStiffnessMatrix();

};

/** u-p Mixed formulation
 * Specialisation for incompressible solid mechanics, mixed formulation without static condensation.
 * The solverSolutionVariable contains the D components of displacements based on HighOrderFunctionSpaceType, then 1 component of pressure based on LowOrderFunctionSpaceType.
 */
template<typename LowOrderFunctionSpaceType, typename HighOrderFunctionSpaceType, typename MixedQuadratureType, typename Term>
class FiniteElementMethodMatrix<
  FunctionSpace::Mixed<LowOrderFunctionSpaceType,HighOrderFunctionSpaceType>,
  MixedQuadratureType,
  Term,
  std::enable_if_t<LowOrderFunctionSpaceType::BasisFunction::isNodalBased, typename HighOrderFunctionSpaceType::Mesh>,
  Equation::isIncompressible<Term>
> :
  public SolidMechanicsCommon<FunctionSpace::Mixed<LowOrderFunctionSpaceType,HighOrderFunctionSpaceType>, HighOrderFunctionSpaceType, MixedQuadratureType, Term>  // this inherits from SolidMechanicsNonlinearSolve and FiniteElementMethodBase<FunctionSpaceType, MixedQuadratureType, Term>
{
public:
  typedef FunctionSpace::Mixed<LowOrderFunctionSpaceType,HighOrderFunctionSpaceType> MixedFunctionSpace;

  // use constructor of base class
  using SolidMechanicsCommon<FunctionSpace::Mixed<LowOrderFunctionSpaceType,HighOrderFunctionSpaceType>, HighOrderFunctionSpaceType, MixedQuadratureType, Term>::SolidMechanicsCommon;

  //! assemble the stiffness matrix
  void setStiffnessMatrix(std::shared_ptr<PartitionedPetscMat<HighOrderFunctionSpaceType>> stiffnessMatrix);

  //! set the internal displacements and pressure variables as copy of the given values
  void setFromSolverVariableSolution(Vec &solverSolutionVariable);

  //! evaluate the nonlinear function that is going to be solved with the Newton scheme
  void evaluateNonlinearFunction(Vec &result);

protected:

  //! get the number of unknowns, also counting displacement values for which Dirichlet BC are set as unknown, for mixed formulation sum of u and p unknowns
  int nLocalUnknowns();

  //! For mixed formulation get the pressure values for the element and store them, to be able to compute the interpolated pressure, by getPressure
  void preparePressureInterpolation(element_no_t elementNo) override;

  //! get the pressure in the current element (set previously by preparePressureInterpolation), interpolated for mixed formulation, by constitutive equation from J for penalty formulation
  double getPressure(double deformationGradientDeterminant, VecD<HighOrderFunctionSpaceType::dim()> xi, double &pressureTilde) override;

  //! compute the D_δp Π_L = integral (J-1)*δp dV, store at position in result after displacements
  void computeIncompressibilityConstraint(Vec &result);

  //! get the index offset for the pressure part in the solverSolutionVariable, this is equal to the number of displacements unknown, either reduced or not, depending on this->data_.computeWithReducedVectors()
  dof_no_t getPressureDofOffset();


  std::array<double,LowOrderFunctionSpaceType::nDofsPerElement()> pressureValuesCurrentElement_;  ///< the pressure values of the current element, to be used for interpolating the pressure. They are set by preparePressureInterpolation and used by getPressure
};

/** Penalty
 * Specialisation for incompressible solid mechanics, not mixed formulation, i.e. penalty formulation,
 */
template<typename FunctionSpaceType, typename QuadratureType, typename Term>
class FiniteElementMethodMatrix<
  FunctionSpaceType,
  QuadratureType,
  Term,
  typename FunctionSpaceType::Mesh,
  Equation::isIncompressible<Term>,
  BasisFunction::isNotMixed<typename FunctionSpaceType::BasisFunction>
> :
  public SolidMechanicsCommon<FunctionSpaceType, FunctionSpaceType, QuadratureType, Term>  // this inherits from SolidMechanicsNonlinearSolve and FiniteElementMethodBase<FunctionSpaceType, QuadratureType, Term>
{
public:
  // use constructor of base class
  using SolidMechanicsCommon<FunctionSpaceType, FunctionSpaceType, QuadratureType, Term>::SolidMechanicsCommon;

  //! assemble the stiffness matrix
  void setStiffnessMatrix(std::shared_ptr<PartitionedPetscMat<FunctionSpaceType>> stiffnessMatrix);

  //! set the internal displacement variable as copy of the given values
  void setFromSolverVariableSolution(Vec &solverSolutionVariable);

  //! evaluate the nonlinear function that is going to be solved with the Newton scheme
  void evaluateNonlinearFunction(Vec &result);

protected:

  //! get the number of unknowns, also counting displacement values for which Dirichlet BC are set as unknown, for mixed formulation sum of u and p unknowns
  int nLocalUnknowns();

  //! For mixed formulation get the pressure values for the element and store them, to be able to compute the interpolated pressure, by getPressure
  void preparePressureInterpolation(element_no_t elementNo) override;

  //! get the pressure in the current element (set previously by preparePressureInterpolation), interpolated for mixed formulation, by constitutive equation from J for penalty formulation
  double getPressure(double deformationGradientDeterminant, VecD<FunctionSpaceType::dim()> xi, double &pressureTilde) override;
};

};  // namespace

#include "spatial_discretization/finite_element_method/solid_mechanics/02_stiffness_matrix_incompressible_mixed.tpp"
#include "spatial_discretization/finite_element_method/solid_mechanics/02_stiffness_matrix_incompressible_mixed_condensation.tpp"
#include "spatial_discretization/finite_element_method/solid_mechanics/02_stiffness_matrix_incompressible_penalty.tpp"
