#include "spatial_discretization/finite_element_method/solid_mechanics/02_stiffness_matrix_incompressible.h"

#include <Python.h>  // has to be the first included header
#include <memory>
#include <petscksp.h>
#include <petscsys.h>

#include "semt/Semt.h"
#include "semt/Shortcuts.h"

#ifdef QUADRATURE_TEST
#include "quadrature/gauss.h"
#endif

namespace SpatialDiscretization
{

// set stiffness matrix for a u-p mixed formulation in which the pressure is not condensed out
template<typename LowOrderBasisOnMeshType, typename HighOrderBasisOnMeshType, typename MixedQuadratureType, typename Term>
void FiniteElementMethodMatrix<
  BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>,
  MixedQuadratureType,
  Term,
  std::enable_if_t<LowOrderBasisOnMeshType::BasisFunction::isNodalBased, typename HighOrderBasisOnMeshType::Mesh>,
  Equation::isIncompressible<Term>
>::
setStiffnessMatrix(Mat stiffnessMatrix)
{
  Mat &tangentStiffnessMatrix = (stiffnessMatrix == PETSC_NULL ? this->data_.tangentStiffnessMatrix() : stiffnessMatrix);

  LOG(TRACE) << "setStiffnessMatrix";

  // set all non-zero entries of the upper left block, this is the same as in penalty formulation and therefore implemented in SolidMechanicsCommon
  this->setStiffnessMatrixEntriesForDisplacements(tangentStiffnessMatrix);

  // set entries for the off-diagonal blocks

  // get pointer to mesh object
  std::shared_ptr<LowOrderBasisOnMeshType> meshP = this->data_.mixedMesh()->lowOrderBasisOnMesh();
  std::shared_ptr<HighOrderBasisOnMeshType> meshU = this->data_.mixedMesh()->highOrderBasisOnMesh();

  const int D = LowOrderBasisOnMeshType::dim();  // = 2 or 3
  const int nElements = meshP->nElements();
  const int nPressureDofsPerElement = LowOrderBasisOnMeshType::nDofsPerElement();
  const int nDisplacementsDofsPerElement = HighOrderBasisOnMeshType::nDofsPerElement();
  const int nDisplacementsUnknownsPerElement = nDisplacementsDofsPerElement * D;

  // define shortcuts for quadrature
  typedef typename MixedQuadratureType::LowOrderQuadrature QuadratureType;
  typedef Quadrature::TensorProduct<D,QuadratureType> QuadratureDD;

  // define type to hold evaluations of integrand
  typedef MathUtility::Matrix<nPressureDofsPerElement,nDisplacementsUnknownsPerElement> EvaluationsType;
  typedef std::array<
            EvaluationsType,
            QuadratureDD::numberEvaluations()
          > EvaluationsArrayType;     // evaluations[nGP^D](nUnknows p,nUnknows u)

  // setup arrays used for integration
  std::array<VecD<D>, QuadratureDD::numberEvaluations()> samplingPoints = QuadratureDD::samplingPoints();
  EvaluationsArrayType evaluationsArray{};

#ifdef QUADRATURE_TEST

  typedef EXACT_QUADRATURE QuadratureExactType;
  typedef Quadrature::TensorProduct<D,QuadratureExactType> QuadratureExactDD;

  // define type to hold evaluations of integrand
  typedef std::array<
            EvaluationsType,
            QuadratureExactDD::numberEvaluations()
          > EvaluationsExactArrayType;     // evaluations[nGP^D](nUnknows p,nUnknows u)

  // setup arrays used for integration
  std::array<VecD<D>, QuadratureExactDD::numberEvaluations()> samplingPointsExact = QuadratureExactDD::samplingPoints();
  EvaluationsExactArrayType evaluationsExactArray{};
#endif

  PetscErrorCode ierr;
  const int pressureDofOffset = this->data_.mixedMesh()->highOrderBasisOnMesh()->nDofs()*D;

  // set values to zero
  // loop over elements
  for (int elementNo = 0; elementNo < nElements; elementNo++)
  {
   // get indices of element-local dofs
    std::array<dof_no_t,nDisplacementsDofsPerElement> dofNoU = meshU->getElementDofNos(elementNo);
    std::array<dof_no_t,nPressureDofsPerElement> dofNoP = meshP->getElementDofNos(elementNo);

    // loop over indices of unknows
    for (int dofIndexL = 0; dofIndexL < nPressureDofsPerElement; dofIndexL++)           // index over pressure dofs in element, L in derivation
    {
      for (int dofIndexM = 0; dofIndexM < nDisplacementsDofsPerElement; dofIndexM++)           // index over displacement dofs in element, M in derivation
      {
        for (int componentIndex = 0; componentIndex < D; componentIndex++)     // index over the component of the dof, a in derivation
        {
          // Compute indices in stiffness matrix
          dof_no_t matrixRowIndex = pressureDofOffset + dofNoP[dofIndexL];
          dof_no_t matrixColumnIndex = dofNoU[dofIndexM]*D + componentIndex;

          VLOG(3) << "set tangentStiffnessMatrix entries " << matrixRowIndex << ", " << matrixColumnIndex;

          // set two symmetric entries
          ierr = MatSetValue(tangentStiffnessMatrix, matrixRowIndex, matrixColumnIndex, 0.0, INSERT_VALUES); CHKERRV(ierr);
          ierr = MatSetValue(tangentStiffnessMatrix, matrixColumnIndex, matrixRowIndex, 0.0, INSERT_VALUES); CHKERRV(ierr);
        }  // a
      }  // M
    }  // L
  }  // elementNo

  // loop over elements
  for (int elementNo = 0; elementNo < nElements; elementNo++)
  {

#ifdef QUADRATURE_TEST
    Control::PerformanceMeasurement::start("stiffnessMatrixPressure");
#endif

    // get geometry field of reference configuration, note the dimension of the vecs is always 3 also for 2D problems
    std::array<Vec3,nDisplacementsDofsPerElement> geometryReferenceValues;
    this->data_.geometryReference().getElementValues(elementNo, geometryReferenceValues);

    // get displacement field values for element
    std::array<VecD<D>,nDisplacementsDofsPerElement> displacementValues;
    this->data_.displacements().getElementValues(elementNo, displacementValues);

    // loop over integration points (e.g. gauss points) for displacement field
    for (unsigned int samplingPointIndex = 0; samplingPointIndex < samplingPoints.size(); samplingPointIndex++)
    {
      // get parameter values of current sampling point
      VecD<D> xi = samplingPoints[samplingPointIndex];

      // compute the 3xD jacobian of the parameter space to world space mapping
      Tensor2<D> jacobianMaterial = MathUtility::transformToDxD<D,D>(HighOrderBasisOnMeshType::computeJacobian(geometryReferenceValues, xi));
      double jacobianDeterminant;
      Tensor2<D> inverseJacobianMaterial = MathUtility::computeInverse<D>(jacobianMaterial, jacobianDeterminant);
      // jacobianMaterial[columnIdx][rowIdx] = dX_rowIdx/dxi_columnIdx
      // inverseJacobianMaterial[columnIdx][rowIdx] = dxi_rowIdx/dX_columnIdx because of inverse function theorem

      // F
      Tensor2<D> deformationGradient = this->computeDeformationGradient(displacementValues, inverseJacobianMaterial, xi);
      double deformationGradientDeterminant = MathUtility::computeDeterminant<D>(deformationGradient);  // J

      VLOG(2) << "";
      VLOG(2) << "element " << elementNo << " xi: " << xi;
      VLOG(2) << "  geometryReferenceValues: " << geometryReferenceValues;
      VLOG(2) << "  displacementValues: " << displacementValues;
      VLOG(2) << "  Jacobian: J_phi=" << jacobianMaterial;
      VLOG(2) << "  deformationGradientDeterminant: J=" << deformationGradientDeterminant;

      std::array<VecD<D>,nDisplacementsDofsPerElement> gradPhi = meshU->getGradPhi(xi);
      // (column-major storage) gradPhi[M][a] = dphi_M / dxi_a
      // gradPhi[column][row] = gradPhi[dofIndex][i] = dphi_dofIndex/dxi_i, columnIdx = dofIndex, rowIdx = which direction

      // loop basis functions and evaluate integrand at xi
      for (int dofIndexL = 0; dofIndexL < nPressureDofsPerElement; dofIndexL++)           // index over pressure dofs in element, L in derivation
      {
        const double psiL = meshP->phi(dofIndexL, xi);

        for (int dofIndexM = 0; dofIndexM < nDisplacementsDofsPerElement; dofIndexM++)           // index over displacement dofs in element, M in derivation
        {
          for (int componentIndex = 0; componentIndex < D; componentIndex++)     // index over the component of the dof, a in derivation
          {
            const double dphiM_da = gradPhi[dofIndexM][componentIndex];   // note that dphi^M_a = dphi^M, i.e. dphi^M_{b,a} = dphi^M_{c,a} = dphi^M_{,a}

            const double integrand = deformationGradientDeterminant * psiL * dphiM_da;

            // store integrand in evaluations array
            evaluationsArray[samplingPointIndex](dofIndexL,dofIndexM*D + componentIndex) = integrand;
          }
        }  // a
      }  // L
    }  // function evaluations

    // integrate all values for result vector entries at once
    EvaluationsType integratedValues = QuadratureDD::computeIntegral(evaluationsArray);

#ifdef QUADRATURE_TEST
    Control::PerformanceMeasurement::stop("stiffnessMatrixPressure");
#endif

#ifdef QUADRATURE_TEST

    // loop over integration points (e.g. gauss points) for displacement field
    for (unsigned int samplingPointIndex = 0; samplingPointIndex < samplingPointsExact.size(); samplingPointIndex++)
    {
      // get parameter values of current sampling point
      VecD<D> xi = samplingPointsExact[samplingPointIndex];

      // compute the 3xD jacobian of the parameter space to world space mapping
      Tensor2<D> jacobianMaterial = MathUtility::transformToDxD<D,D>(HighOrderBasisOnMeshType::computeJacobian(geometryReferenceValues, xi));
      double jacobianDeterminant;
      Tensor2<D> inverseJacobianMaterial = MathUtility::computeInverse<D>(jacobianMaterial, jacobianDeterminant);
      // jacobianMaterial[columnIdx][rowIdx] = dX_rowIdx/dxi_columnIdx
      // inverseJacobianMaterial[columnIdx][rowIdx] = dxi_rowIdx/dX_columnIdx because of inverse function theorem

      // F
      Tensor2<D> deformationGradient = this->computeDeformationGradient(displacementValues, inverseJacobianMaterial, xi);
      double deformationGradientDeterminant = MathUtility::computeDeterminant<D>(deformationGradient);  // J

      std::array<VecD<D>,nDisplacementsDofsPerElement> gradPhi = meshU->getGradPhi(xi);
      // (column-major storage) gradPhi[M][a] = dphi_M / dxi_a
      // gradPhi[column][row] = gradPhi[dofIndex][i] = dphi_dofIndex/dxi_i, columnIdx = dofIndex, rowIdx = which direction

      // loop basis functions and evaluate integrand at xi
      for (int dofIndexL = 0; dofIndexL < nPressureDofsPerElement; dofIndexL++)           // index over pressure dofs in element, L in derivation
      {
        const double psiL = meshP->phi(dofIndexL, xi);

        for (int dofIndexM = 0; dofIndexM < nDisplacementsDofsPerElement; dofIndexM++)           // index over displacement dofs in element, M in derivation
        {
          for (int componentIndex = 0; componentIndex < D; componentIndex++)     // index over the component of the dof, a in derivation
          {
            const double dphiM_da = gradPhi[dofIndexM][componentIndex];   // note that dphi^M_a = dphi^M, i.e. dphi^M_{b,a} = dphi^M_{c,a} = dphi^M_{,a}

            const double integrand = deformationGradientDeterminant * psiL * dphiM_da;

            // store integrand in evaluations array
            evaluationsExactArray[samplingPointIndex](dofIndexL,dofIndexM*D + componentIndex) = integrand;
          }
        }  // a
      }  // L
    }  // function evaluations

    // integrate all values for result vector entries at once
    EvaluationsType integratedValuesExact = QuadratureExactDD::computeIntegral(evaluationsExactArray);

    Control::PerformanceMeasurement::measureError("stiffnessMatrixPressure", (integratedValues - integratedValuesExact)/integratedValuesExact);

#endif

    // get indices of element-local dofs
    std::array<dof_no_t,nDisplacementsDofsPerElement> dofNoU = meshU->getElementDofNos(elementNo);
    std::array<dof_no_t,nPressureDofsPerElement> dofNoP = meshP->getElementDofNos(elementNo);

    // add entries in tangent stiffness matrix
    // loop over indices of unknows
    for (int dofIndexL = 0; dofIndexL < nPressureDofsPerElement; dofIndexL++)           // index over pressure dofs in element, L in derivation
    {
      for (int dofIndexM = 0; dofIndexM < nDisplacementsDofsPerElement; dofIndexM++)           // index over displacement dofs in element, M in derivation
      {
        for (int componentIndex = 0; componentIndex < D; componentIndex++)     // index over the component of the dof, a in derivation
        {
          // get integrated value for current dof
          double integratedValue = integratedValues(dofIndexL,dofIndexM*D + componentIndex);

          // Compute indices in stiffness matrix
          dof_no_t matrixRowIndex = pressureDofOffset + dofNoP[dofIndexL];
          dof_no_t matrixColumnIndex = dofNoU[dofIndexM]*D + componentIndex;

          ierr = MatSetValue(tangentStiffnessMatrix, matrixRowIndex, matrixColumnIndex, integratedValue, ADD_VALUES); CHKERRV(ierr);

          // set symmetric entry
          ierr = MatSetValue(tangentStiffnessMatrix, matrixColumnIndex, matrixRowIndex, integratedValue, ADD_VALUES); CHKERRV(ierr);

        }  // a
      }  // M
    }  // L
  }  // elementNo

  // because this is used in nonlinear solver context, assembly has to be called here, not via data->finalAssembly
  ierr = MatAssemblyBegin(tangentStiffnessMatrix, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);
  ierr = MatAssemblyEnd(tangentStiffnessMatrix, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);

  if (!this->tangentStiffnessMatrixInitialized_)
  {
    VLOG(3) << "tangent stiffness matrix before zero rows columns: " << PetscUtility::getStringMatrix(tangentStiffnessMatrix);
    VLOG(3) << "number dirichletIndices: " << this->dirichletIndices_.size();

    // zero rows and columns for which Dirichlet BC is set
    ierr = MatZeroRowsColumns(tangentStiffnessMatrix, this->dirichletIndices_.size(), this->dirichletIndices_.data(), 1.0, PETSC_IGNORE, PETSC_IGNORE); CHKERRV(ierr);

    VLOG(3) << "tangent stiffness matrix after zero rows columns: " << PetscUtility::getStringMatrix(tangentStiffnessMatrix);

    // set option that all insert/add operations to new nonzero locations will be discarded. This keeps the nonzero structure forever.
    // (The diagonal entries will be set to different values, but that doesn't matter because the Dirichlet values for updates are 0 and thus independent of the diagonal scaling (d*Δx=0 -> Δx=0 independent of d))
    ierr = MatSetOption(tangentStiffnessMatrix, MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE); CHKERRV(ierr);

    this->tangentStiffnessMatrixInitialized_ = true;
  }


  VLOG(1) << "computed analytical stiffness matrix (not reduced): " << PetscUtility::getStringMatrix(tangentStiffnessMatrix);
}

template<typename LowOrderBasisOnMeshType, typename HighOrderBasisOnMeshType, typename MixedQuadratureType, typename Term>
const dof_no_t FiniteElementMethodMatrix<
  BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>,
  MixedQuadratureType,
  Term,
  std::enable_if_t<LowOrderBasisOnMeshType::BasisFunction::isNodalBased, typename HighOrderBasisOnMeshType::Mesh>,
  Equation::isIncompressible<Term>
>::
getPressureDofOffset()
{
  const int D = HighOrderBasisOnMeshType::dim();

  if (this->data_.computeWithReducedVectors())
  {
    return this->data_.mixedMesh()->highOrderBasisOnMesh()->nDofs()*D - this->dirichletIndices_.size();
  }
  else
  {
    return this->data_.mixedMesh()->highOrderBasisOnMesh()->nDofs()*D;
  }
}

template<typename LowOrderBasisOnMeshType, typename HighOrderBasisOnMeshType, typename MixedQuadratureType, typename Term>
void FiniteElementMethodMatrix<
  BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>,
  MixedQuadratureType,
  Term,
  std::enable_if_t<LowOrderBasisOnMeshType::BasisFunction::isNodalBased, typename HighOrderBasisOnMeshType::Mesh>,
  Equation::isIncompressible<Term>
>::
setFromSolverVariableSolution(Vec &solverVariableSolution)
{
  // store displacement values
  const int D = HighOrderBasisOnMeshType::dim();
  const int nDisplacementsUnknowns = this->data_.mixedMesh()->highOrderBasisOnMesh()->nDofs() * D;

  // get memory read access to solverVariableSolution
  const double *solverVariableSolutionData;
  VecGetArrayRead(solverVariableSolution, &solverVariableSolutionData);

  if (this->data_.computeWithReducedVectors())
  {
    // store displacements values
    const int nUnknownsOutputVector = nDisplacementsUnknowns;
    this->expandVector(solverVariableSolution, this->data_.displacements().values(), nUnknownsOutputVector);
  }
  else
  {
    // store displacements values
    double *displacementsData;
    VecGetArray(this->data_.displacements().values(), &displacementsData);

    // copy values from solverVariableSolution to displacements Vec
    for (dof_no_t currentDofNo = 0; currentDofNo < nDisplacementsUnknowns; currentDofNo++)
    {
      displacementsData[currentDofNo] = solverVariableSolutionData[currentDofNo];
    }

    // return memory access to Petsc
    VecRestoreArray(this->data_.displacements().values(), &displacementsData);
  }

  // store pressure values
  double *pressureData;
  VecGetArray(this->data_.pressure().values(), &pressureData);

  const int nDofsPressure = this->data_.mixedMesh()->lowOrderBasisOnMesh()->nDofs();
  const int pressureDofOffset = getPressureDofOffset();

  for (dof_no_t currentDofNo = 0; currentDofNo < nDofsPressure; currentDofNo++)
  {
    pressureData[currentDofNo] = solverVariableSolutionData[pressureDofOffset + currentDofNo];
  }
  VecRestoreArray(this->data_.pressure().values(), &pressureData);

  VecRestoreArrayRead(solverVariableSolution, &solverVariableSolutionData);

  if (VLOG_IS_ON(1))
  {
    VLOG(1) << "extracted displacements: " << PetscUtility::getStringVector(this->data_.displacements().values());
    VLOG(1) << "extracted pressure: " << PetscUtility::getStringVector(this->data_.pressure().values());
  }
}

template<typename LowOrderBasisOnMeshType, typename HighOrderBasisOnMeshType, typename MixedQuadratureType, typename Term>
void FiniteElementMethodMatrix<
  BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>,
  MixedQuadratureType,
  Term,
  std::enable_if_t<LowOrderBasisOnMeshType::BasisFunction::isNodalBased, typename HighOrderBasisOnMeshType::Mesh>,
  Equation::isIncompressible<Term>
>::
evaluateNonlinearFunction(Vec &result)
{
  this->computeInternalMinusExternalVirtualWork(result);
  this->computeIncompressibilityConstraint(result);
}

template<typename LowOrderBasisOnMeshType, typename HighOrderBasisOnMeshType, typename MixedQuadratureType, typename Term>
void FiniteElementMethodMatrix<
  BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>,
  MixedQuadratureType,
  Term,
  std::enable_if_t<LowOrderBasisOnMeshType::BasisFunction::isNodalBased, typename HighOrderBasisOnMeshType::Mesh>,
  Equation::isIncompressible<Term>
>::
computeIncompressibilityConstraint(Vec &result)
{
  // compute int_Omega (J(u) - 1) δp dV
  LOG(TRACE) << "computeIncompressibilityConstraint";

  // get pointer to mesh object
  std::shared_ptr<LowOrderBasisOnMeshType> mesh = this->data_.mixedMesh()->lowOrderBasisOnMesh();

  const int D = LowOrderBasisOnMeshType::dim();  // = 2 or 3
  const int nElements = mesh->nElements();
  const int nPressureDofsPerElement = LowOrderBasisOnMeshType::nDofsPerElement();
  const int nDisplacementsDofsPerElement = HighOrderBasisOnMeshType::nDofsPerElement();

  // define shortcuts for quadrature
  typedef typename MixedQuadratureType::LowOrderQuadrature QuadratureType;
  typedef Quadrature::TensorProduct<D,QuadratureType> QuadratureDD;

  // define type to hold evaluations of integrand
  typedef std::array<double,nPressureDofsPerElement> EvaluationsType;
  typedef std::array<
            EvaluationsType,
            QuadratureDD::numberEvaluations()
          > EvaluationsArrayType;     // evaluations[nGP^D](nUnknows,nUnknows)

  // setup arrays used for integration
  std::array<VecD<D>, QuadratureDD::numberEvaluations()> samplingPoints = QuadratureDD::samplingPoints();
  EvaluationsArrayType evaluationsArray{};

  // get memory from PETSc where to store result
  double *resultData;
  VecGetArray(result, &resultData);

  // get correct offset
  const dof_no_t pressureDofOffset = getPressureDofOffset();

  // set all values corresponding to the incompressiblityt constraint to zero
  double *vectorBegin = pressureDofOffset + resultData;
  double *vectorEnd = vectorBegin + this->data_.mixedMesh()->lowOrderBasisOnMesh()->nDofs();
  std::fill(vectorBegin, vectorEnd, 0.0);

  // loop over elements
  for (int elementNo = 0; elementNo < nElements; elementNo++)
  {
    // get geometry field of reference configuration, note the dimension of the vecs is always 3 also for 2D problems
    std::array<Vec3,nDisplacementsDofsPerElement> geometryReferenceValues;
    this->data_.geometryReference().getElementValues(elementNo, geometryReferenceValues);

    // get displacement field values for element
    std::array<VecD<D>,nDisplacementsDofsPerElement> displacementValues;
    this->data_.displacements().getElementValues(elementNo, displacementValues);

    // loop over integration points (e.g. gauss points) for displacement field
    for (unsigned int samplingPointIndex = 0; samplingPointIndex < samplingPoints.size(); samplingPointIndex++)
    {
      // get parameter values of current sampling point
      VecD<D> xi = samplingPoints[samplingPointIndex];

      // compute the 3xD jacobian of the parameter space to world space mapping
      Tensor2<D> jacobianMaterial = MathUtility::transformToDxD<D,D>(HighOrderBasisOnMeshType::computeJacobian(geometryReferenceValues, xi));
      double jacobianDeterminant;
      Tensor2<D> inverseJacobianMaterial = MathUtility::computeInverse<D>(jacobianMaterial, jacobianDeterminant);
      // jacobianMaterial[columnIdx][rowIdx] = dX_rowIdx/dxi_columnIdx
      // inverseJacobianMaterial[columnIdx][rowIdx] = dxi_rowIdx/dX_columnIdx because of inverse function theorem

      checkInverseIsCorrect<D>(jacobianMaterial, inverseJacobianMaterial, "jacobianMaterial");

      // F
      Tensor2<D> deformationGradient = this->computeDeformationGradient(displacementValues, inverseJacobianMaterial, xi);
      double deformationGradientDeterminant = MathUtility::computeDeterminant<D>(deformationGradient);  // J

      VLOG(2) << "";
      VLOG(2) << "element " << elementNo << " xi: " << xi;
      VLOG(2) << "  geometryReferenceValues: " << geometryReferenceValues;
      VLOG(2) << "  displacementValues: " << displacementValues;
      VLOG(2) << "  Jacobian: J_phi=" << jacobianMaterial;
      VLOG(2) << "  deformationGradientDeterminant: J=" << deformationGradientDeterminant;

      // loop basis functions and evaluate integrand at xi
      for (int dofIndex = 0; dofIndex < nPressureDofsPerElement; dofIndex++)           // index over dofs in element, L in derivation
      {
        const double phiL = mesh->phi(dofIndex, xi);
        const double integrand = (deformationGradientDeterminant - 1.0) * phiL;

        // store integrand in evaluations array
        evaluationsArray[samplingPointIndex][dofIndex] = integrand;
      }  // L
    }  // function evaluations

    // integrate all values for result vector entries at once
    EvaluationsType integratedValues = QuadratureDD::computeIntegral(evaluationsArray);

    // get indices of element-local dofs
    std::array<dof_no_t,nPressureDofsPerElement> dofNo = mesh->getElementDofNos(elementNo);

    // add entries in result vector
    // loop over indices of unknows
    for (int dofIndex = 0; dofIndex < nPressureDofsPerElement; dofIndex++)           // index over dofs in element, L in derivation
    {
      // get integrated value for current dof
      double integratedValue = integratedValues[dofIndex];

      // compute index in result vector
      dof_no_t resultVectorIndex = dofNo[dofIndex];

      // store value to result Vec
      resultData[pressureDofOffset + resultVectorIndex] += integratedValue;

    }  // dofIndex
  }  // elementNo

  if (VLOG_IS_ON(1))
  {
    std::stringstream s;
    s << "(" << resultData[pressureDofOffset + 0];

    for (int i = 1; i < this->data_.mixedMesh()->lowOrderBasisOnMesh()->nDofs(); i++)
    {
      s << "," << resultData[pressureDofOffset + i];
    }
    s << ")";
    VLOG(1) << "  incomp: " << s.str();
  }

  // return memory access of result vector back to PETSc
  VecRestoreArray(result, &resultData);

}

template<typename LowOrderBasisOnMeshType, typename HighOrderBasisOnMeshType, typename MixedQuadratureType, typename Term>
void FiniteElementMethodMatrix<
  BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>,
  MixedQuadratureType,
  Term,
  std::enable_if_t<LowOrderBasisOnMeshType::BasisFunction::isNodalBased, typename HighOrderBasisOnMeshType::Mesh>,
  Equation::isIncompressible<Term>
>::
preparePressureInterpolation(element_no_t elementNo)
{
  // get pressure field values for element
  this->data_.pressure().getElementValues(elementNo, pressureValuesCurrentElement_);
}

template<typename LowOrderBasisOnMeshType, typename HighOrderBasisOnMeshType, typename MixedQuadratureType, typename Term>
double FiniteElementMethodMatrix<
  BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>,
  MixedQuadratureType,
  Term,
  std::enable_if_t<LowOrderBasisOnMeshType::BasisFunction::isNodalBased, typename HighOrderBasisOnMeshType::Mesh>,
  Equation::isIncompressible<Term>
>::
getPressure(double deformationGradientDeterminant, VecD<HighOrderBasisOnMeshType::dim()> xi, double &pressureTilde)
{
  // interpolate pressure value at xi in current element
  pressureTilde = this->data_.mixedMesh()->lowOrderBasisOnMesh()->interpolateValueInElement(pressureValuesCurrentElement_, xi);

  VLOG(3) << "pressureValuesCurrentElement: " << pressureValuesCurrentElement_ << ", pressureTilde: " << pressureTilde;

  return pressureTilde;
}

template<typename LowOrderBasisOnMeshType, typename HighOrderBasisOnMeshType, typename MixedQuadratureType, typename Term>
const int FiniteElementMethodMatrix<
  BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>,
  MixedQuadratureType,
  Term,
  std::enable_if_t<LowOrderBasisOnMeshType::BasisFunction::isNodalBased, typename HighOrderBasisOnMeshType::Mesh>,
  Equation::isIncompressible<Term>
>::
nUnknowns()
{
  const int D = HighOrderBasisOnMeshType::dim();
  const int nUnknownsDisplacements = this->data_.mixedMesh()->highOrderBasisOnMesh()->nDofs() * D;
  const int nUnknownsPressure = this->data_.mixedMesh()->lowOrderBasisOnMesh()->nDofs();
  return nUnknownsDisplacements + nUnknownsPressure;
}

};    // namespace