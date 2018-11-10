#pragma once

#include <vector>

#include "spatial_discretization/finite_element_method/solid_mechanics/solid_mechanics_nonlinear_solve.h"
#include "spatial_discretization/finite_element_method/solid_mechanics/solid_mechanics_utility.h"

namespace SpatialDiscretization
{

/** Helper class that provides common methods to the mixed, mixed_condensation and penalty classes.
 *  The FunctionSpaceType is Mixed for mixed formulation, FunctionSpaceTypeForUtility is the class used for SolidMechanicsUtility, which is not a mixed FunctionSpace.
 */
template<typename FunctionSpaceType, typename FunctionSpaceTypeForUtility, typename QuadratureType, typename Term>
class SolidMechanicsCommon :
  public SolidMechanicsNonlinearSolve<FunctionSpaceType, QuadratureType, Term>,  // this inherits from FiniteElementMethodBase<FunctionSpaceType, QuadratureType, Term>
  public SolidMechanicsUtility<FunctionSpaceTypeForUtility, Term>
{
public:

  //! inherited constructor
  using SolidMechanicsNonlinearSolve<FunctionSpaceType, QuadratureType, Term>::SolidMechanicsNonlinearSolve;

  //! compute the tangent stiffness matrix and store it internally
  void setStiffnessMatrix();

  //! return the tangent stiffness matrix, called from a PETSc SNES callback
  std::shared_ptr<PartitionedPetscMat<FunctionSpaceType>> tangentStiffnessMatrix();

  //! return the virtual external energy which is constant
  Vec &rightHandSide();

  //! return the displacements values
  Vec &displacements();

  //! compute the value of dW_int - dW_ext
  void computeInternalMinusExternalVirtualWork(Vec &result);

  //! compute and return the appropriate analytical stiffness matrix
  void computeAnalyticStiffnessMatrix(std::shared_ptr<PartitionedPetscMat<FunctionSpaceType>> solverStiffnessMatrix);

  //! write the current state as output
  void writeOutput();

protected:

  //! compute the tangent stiffnes matrix into the given matrix (implemented by inherited class)
  virtual void setStiffnessMatrix(std::shared_ptr<PartitionedPetscMat<FunctionSpaceType>> stiffnessMatrix) = 0;

  //! set all entries in the given stiffness matrix for the analytic jacobian only with displacements. This corresponds to the full matrix in penalty formulation and to the upper left block in mixed formulation.
  //! do not call final assembly or anything
  void setStiffnessMatrixEntriesForDisplacements(std::shared_ptr<PartitionedPetscMat<FunctionSpaceType>> stiffnessMatrix);

  //! compute the extern virtual work term, dW_ext.
  void computeExternalVirtualWork(Vec &result);

  //! compute the internal virtual work term, dW_int. It is nearly the same for penalty and mixed formulation, the difference is in the pressure compution
  //! (Artifical pressure by constitutive assumption in penalty formulation, separately interpolated pressure in mixed formulation) This difference is realized by the callback function to getPressure(xi)
  void computeInternalVirtualWork(Vec &resultVec);

  //! get the extern virtual work term, dW_ext. It is computed if it depends on displacements, otherwise the stored value is returned.
  void getExternalVirtualWork(Vec &result);

  //! initialize everything, set rhs and compute tangent stiffness matrix
  virtual void initialize();

  //! read material parameters from config and set the values for static expressions within SEMT
  void initializeMaterialParameters();

  //! return the bool value of data_, whether to do the nonlinear solution completely in reduced vectors without the Dirichlet BC variables
  bool computeWithReducedVectors();

  //! For mixed formulation get the pressure values for the element and store them, to be able to compute the interpolated pressure, by getPressure
  virtual void preparePressureInterpolation(element_no_t elementNo) = 0;

  //! get the pressure in the current element (set previously by preparePressureInterpolation), interpolated for mixed formulation, by constitutive equation from J for penalty formulation
  virtual double getPressure(double deformationGradientDeterminant, VecD<FunctionSpaceType::dim()> xi, double &pressureTilde) = 0;

  bool tangentStiffnessMatrixInitialized_ = false;   ///< if the tangentStiffnessMatrix has been set
  bool outputIntermediateSteps_ = false;    ///< if intermediate steps while solving should be output
};

};  // namespace

#include "spatial_discretization/finite_element_method/solid_mechanics/solid_mechanics_common.tpp"
