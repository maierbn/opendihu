#include "specialized_solver/hyperelasticity/quasi_static_hyperelasticity_solver.h"

#include <Python.h>  // has to be the first included header

#include "equation/mooney_rivlin_incompressible.h"

namespace TimeSteppingScheme
{

void QuasiStaticHyperelasticitySolver::
materialComputeResidual()
{
  // compute Wint - Wext in solverVariableResidual_
  // if useNestedMat_:
  //  output is solverVariableResidual_, which is a nested vector of subvectorsResidual_[0-4]
  //  input is solverVariableSolution_, which is a nested vector of subvectorsSolution_[0-4]
  //  subvectorsSolution_ is the same vector as this->data_.displacements() (3 components) and this->data_.pressure() (1 component)
  //
  // if not useNestedMat_:
  //  output is solverVariableResidual_, a normal Vec, no Dirichlet BC dofs
  //  input is solverVariableSolution_, a normal Vec, the same values have already been assigned to this->data_.displacements() and this->data_.pressure()

  // -----------------------------
  // compute internal virtual work
  // δW_int = int_Ω 1/2 S_AB (F_aB phi_L,A + F_aA phi_L,B) dV

  //LOG(TRACE) << "materialComputeResidual";
  const bool outputValues = false;
  if (outputValues)
    LOG(DEBUG) << "input: " << getString(solverVariableSolution_);

  assert(combinedVecResidual_->currentRepresentation() == Partition::values_representation_t::representationCombinedGlobal);
  assert(combinedVecSolution_->currentRepresentation() == Partition::values_representation_t::representationCombinedGlobal);

  static int evaluationNo = 0;  // counter how often this function was called

  // get pointer to function space
  std::shared_ptr<DisplacementsFunctionSpace> displacementsFunctionSpace = this->data_.displacementsFunctionSpace();
  std::shared_ptr<PressureFunctionSpace> pressureFunctionSpace = this->data_.pressureFunctionSpace();

  const int D = 3;  // dimension
  const int nDisplacementsDofsPerElement = DisplacementsFunctionSpace::nDofsPerElement();
  const int nPressureDofsPerElement = PressureFunctionSpace::nDofsPerElement();
  const int nElementsLocal = displacementsFunctionSpace->nElementsLocal();
  const int nUnknowsPerElement = nDisplacementsDofsPerElement*D;    // D directions for displacements per dof

  // define shortcuts for quadrature
  typedef Quadrature::TensorProduct<D,Quadrature::Gauss<3>> QuadratureDD;   // quadratic*quadratic = 4th order polynomial, 3 gauss points = 2*3-1 = 5th order exact

  // define type to hold evaluations of integrand
  typedef std::array<double, nUnknowsPerElement> EvaluationsDisplacementsType;
  std::array<EvaluationsDisplacementsType, QuadratureDD::numberEvaluations()> evaluationsArrayDisplacements{};

  typedef std::array<double, nPressureDofsPerElement> EvaluationsPressureType;
  std::array<EvaluationsPressureType, QuadratureDD::numberEvaluations()> evaluationsArrayPressure{};

  // setup arrays used for integration
  std::array<Vec3, QuadratureDD::numberEvaluations()> samplingPoints = QuadratureDD::samplingPoints();

  PetscErrorCode ierr;

  // set values to zero
  if (useNestedMat_)
  {
    ierr = VecZeroEntries(solverVariableResidual_); CHKERRV(ierr);
  }
  else
  {
    combinedVecResidual_->zeroEntries();
    combinedVecResidual_->startGhostManipulation();
  }

  // dump input vector
  // dumpVector(std::string filename, std::string format, Vec &vector, MPI_Comm mpiCommunicator, int componentNo=0, int nComponents=1);
  std::stringstream filename;
  filename << "out/x" << std::setw(3) << std::setfill('0') << evaluationNo;
  //PetscUtility::dumpVector(filename.str(), "matlab", solverVariableSolution_, displacementsFunctionSpace->meshPartition()->mpiCommunicator());
  combinedVecSolution_->dumpGlobalNatural(filename.str());

  // loop over elements
  for (int elementNoLocal = 0; elementNoLocal < nElementsLocal; elementNoLocal++)
  {
    // get geometry field of reference configuration
    std::array<Vec3,nDisplacementsDofsPerElement> geometryReferenceValues;
    this->data_.geometryReference()->getElementValues(elementNoLocal, geometryReferenceValues);

    // get displacements field values for element
    std::array<Vec3,nDisplacementsDofsPerElement> displacementsValues;
    this->data_.displacements()->getElementValues(elementNoLocal, displacementsValues);

    std::array<double,nPressureDofsPerElement> pressureValuesCurrentElement;
    this->data_.pressure()->getElementValues(elementNoLocal, pressureValuesCurrentElement);

    // loop over integration points (e.g. gauss points) for displacements field
    for (unsigned int samplingPointIndex = 0; samplingPointIndex < samplingPoints.size(); samplingPointIndex++)
    {
      // get parameter values of current sampling point
      Vec3 xi = samplingPoints[samplingPointIndex];

      // compute the 3x3 jacobian of the parameter space to world space mapping
      Tensor2<D> jacobianMaterial = DisplacementsFunctionSpace::computeJacobian(geometryReferenceValues, xi);
      double jacobianDeterminant;
      Tensor2<D> inverseJacobianMaterial = MathUtility::computeInverse<D>(jacobianMaterial, jacobianDeterminant);

      // jacobianMaterial[columnIdx][rowIdx] = dX_rowIdx/dxi_columnIdx
      // inverseJacobianMaterial[columnIdx][rowIdx] = dxi_rowIdx/dX_columnIdx because of inverse function theorem

      // F
      Tensor2<D> deformationGradient = this->computeDeformationGradient(displacementsValues, inverseJacobianMaterial, xi);
      double deformationGradientDeterminant = MathUtility::computeDeterminant<D>(deformationGradient);  // J

      Tensor2<D> rightCauchyGreen = this->computeRightCauchyGreenTensor(deformationGradient);  // C = F^T*F

      double rightCauchyGreenDeterminant;   // J^2
      Tensor2<D> inverseRightCauchyGreen = MathUtility::computeSymmetricInverse<D>(rightCauchyGreen, rightCauchyGreenDeterminant);  // C^-1

      // invariants
      std::array<double,3> invariants = this->computeInvariants(rightCauchyGreen, rightCauchyGreenDeterminant);  // I_1, I_2, I_3
      std::array<double,2> reducedInvariants = this->computeReducedInvariants(invariants, deformationGradientDeterminant); // Ibar_1, Ibar_2

      // pressure is the separately interpolated pressure for mixed formulation
      double pressure = pressureFunctionSpace->interpolateValueInElement(pressureValuesCurrentElement, xi);

      // Pk2 stress tensor S = S_vol + S_iso (p.234)
      //! compute 2nd Piola-Kirchhoff stress tensor S = 2*dPsi/dC and the fictitious PK2 Stress Sbar
      Tensor2<D> PK2Stress = this->computePK2Stress(pressure, rightCauchyGreen, inverseRightCauchyGreen, reducedInvariants, deformationGradientDeterminant);

      std::array<Vec3,nDisplacementsDofsPerElement> gradPhi = displacementsFunctionSpace->getGradPhi(xi);
      // (column-major storage) gradPhi[L][a] = dphi_L / dxi_a
      // gradPhi[column][row] = gradPhi[dofIndex][i] = dphi_dofIndex/dxi_i, columnIdx = dofIndex, rowIdx = which direction

      VLOG(2) << "";
      VLOG(2) << "element " << elementNoLocal << " xi: " << xi;
      VLOG(2) << "  geometryReferenceValues: " << geometryReferenceValues;
      VLOG(2) << "  displacementsValues: " << displacementsValues;
      VLOG(2) << "  Jacobian: J_phi=" << jacobianMaterial;
      VLOG(2) << "  jacobianDeterminant: J=" << jacobianDeterminant;
      VLOG(2) << "  inverseJacobianMaterial: J_phi^-1=" << inverseJacobianMaterial;
      VLOG(2) << "  deformationGradient: F=" << deformationGradient;
      VLOG(2) << "  deformationGradientDeterminant: det F=" << deformationGradientDeterminant;
      VLOG(2) << "  rightCauchyGreen: C=" << rightCauchyGreen;
      VLOG(2) << "  rightCauchyGreenDeterminant: det C=" << rightCauchyGreenDeterminant;
      VLOG(2) << "  inverseRightCauchyGreen: C^-1=" << inverseRightCauchyGreen;
      VLOG(2) << "  invariants: I1,I2,I3: " << invariants;
      VLOG(2) << "  reducedInvariants: Ibar1, Ibar2: " << reducedInvariants;
      VLOG(2) << "  pressure/artificialPressure: " << pressure;
      //VLOG(2) << "  artificialPressure: p=" << artificialPressure << ", artificialPressureTilde: pTilde=" << artificialPressureTilde;
      VLOG(2) << "  PK2Stress: S=" << PK2Stress;
      VLOG(2) << "  gradPhi: " << gradPhi;

      VLOG(1) << "  xi: " << xi << ", J: " << deformationGradientDeterminant << ", p: " << pressure << ", S11: " << PK2Stress[0][0];

      if (samplingPointIndex == 0 && D == 3)
        VLOG(1) << " F11: " << deformationGradient[0][0] << ", F22,F33: " << deformationGradient[1][1] << "," << deformationGradient[2][2]
          << ", F12,F13,F23: " << deformationGradient[1][0] << "," << deformationGradient[2][0] << "," << deformationGradient[2][1]
          << ", J: " << deformationGradientDeterminant << ", p: " << pressure;

      /*if (VLOG_IS_ON(2))
      {
        Tensor2<D> greenLangrangeStrain = this->computeGreenLagrangeStrain(rightCauchyGreen);
        VLOG(2) << "  strain E=" << greenLangrangeStrain;
      }*/

      if (deformationGradientDeterminant < 1e-12)   // if any entry of the deformation gradient is negative
      {
        LOG(FATAL) << "Deformation gradient " << deformationGradient << " has zero or negative determinant " << deformationGradientDeterminant
          << std::endl << "Geometry values in element " << elementNoLocal << ": " << geometryReferenceValues << std::endl
          << "Displacements at xi " << xi << ": " << displacementsValues;
      }

      // loop over basis functions and evaluate integrand at xi for displacement part (δW_int - δW_ext)
      for (int aDof = 0; aDof < nDisplacementsDofsPerElement; aDof++)           // index over dofs, each dof has D components, L in derivation
      {
        for (int aComponent = 0; aComponent < D; aComponent++)     // lower-case a in derivation, index over displacements components
        {
          // compute index of degree of freedom and component (result vector index)
          const int i = D*aDof + aComponent;

          // compute result[i]
          double integrand = 0.0;
          for (int aInternal = 0; aInternal < D; aInternal++)     // capital A in derivation
          {
            for (int bInternal = 0; bInternal < D; bInternal++)     // capital B in derivation
            {
              const double faB = deformationGradient[bInternal][aComponent];
              const double faA = deformationGradient[aInternal][aComponent];

              // ----------------------------
              // compute derivatives of phi
              double dphiL_dXA = 0.0;
              double dphiL_dXB = 0.0;

              // helper index k for multiplication with inverse Jacobian
              for (int k = 0; k < D; k++)
              {
                // compute dphiL/dXA from dphiL/dxik and dxik/dXA
                const double dphiL_dxik = gradPhi[aDof][k];    // dphi_L/dxik
                const double dxik_dXA = inverseJacobianMaterial[aInternal][k];  // inverseJacobianMaterial[A][k] = J^{-1}_kA = dxi_k/dX_A

                dphiL_dXA += dphiL_dxik * dxik_dXA;

                // compute dphiL/dXB from dphiL/dxik and dxik/dXB
                const double dxik_dXB = inverseJacobianMaterial[bInternal][k];  // inverseJacobianMaterial[B][k] = J^{-1}_kB = dxi_k/dX_B

                dphiL_dXB += dphiL_dxik * dxik_dXB;
              }

              integrand += 1./2. * PK2Stress[bInternal][aInternal] * (faB*dphiL_dXA + faA*dphiL_dXB);

            }  // B, bInternal
          }  // A, aInternal

          VLOG(2) << "   (L,a)=(" << aDof << "," << aComponent << "), integrand: " << integrand;

          // store integrand in evaluations array
          evaluationsArrayDisplacements[samplingPointIndex][i] = integrand;

        }  // a
      }  // L

      // loop over basis functions and evaluate integrand at xi for pressure part ((J-1)*psi)
      for (int dofIndex = 0; dofIndex < nPressureDofsPerElement; dofIndex++)           // index over dofs in element, L in derivation
      {
        const double phiL = pressureFunctionSpace->phi(dofIndex, xi);
        const double integrand = (deformationGradientDeterminant - 1.0) * phiL;     // (J-1) * phi_L

        // store integrand in evaluations array
        evaluationsArrayPressure[samplingPointIndex][dofIndex] = integrand;
      }  // L

    }  // function evaluations

    // integrate all values for result vector entries at once
    EvaluationsDisplacementsType integratedValuesDisplacements = QuadratureDD::computeIntegral(evaluationsArrayDisplacements);
    EvaluationsPressureType integratedValuesPressure = QuadratureDD::computeIntegral(evaluationsArrayPressure);

    // get indices of element-local dofs
    std::array<dof_no_t,nDisplacementsDofsPerElement> dofNosLocal = displacementsFunctionSpace->getElementDofNosLocal(elementNoLocal);
    std::array<dof_no_t,nPressureDofsPerElement> pressureDofNosLocal = pressureFunctionSpace->getElementDofNosLocal(elementNoLocal);

    VLOG(2) << "  element " << elementNoLocal << " has dofs " << dofNosLocal;

    // add entries in result vector for displacements
    // loop over indices of unknows (aDof,aComponent)
    for (int aDof = 0; aDof < nDisplacementsDofsPerElement; aDof++)      // L, this is the dof within the element
    {
      for (int aComponent = 0; aComponent < D; aComponent++)    // a
      {
        // compute index of degree of freedom and component (matrix row index)
        const int i = D*aDof + aComponent;

        // integrate value and set entry
        double integratedValue = integratedValuesDisplacements[i];

        // get local dof no, aDof is the dof within the element, dofNoLocal is the dof within the local subdomain
        dof_no_t dofNoLocal = dofNosLocal[aDof];

        VLOG(1) << "  result vector (L,a)=(" <<aDof<< "," <<aComponent<< "), " <<i<< ", dof " << dofNosLocal[aDof] << " all elemental dofs: "  << dofNosLocal
          << ", integrated value: " <<integratedValue;

        if (useNestedMat_)
        {
          ierr = VecSetValue(subvectorsResidual_[aComponent], dofNoLocal, integratedValue, ADD_VALUES); CHKERRV(ierr);
        }
        else
        {
          combinedVecResidual_->setValue(aComponent, dofNoLocal, integratedValue, ADD_VALUES);
          VLOG(1) << "u: set value " << integratedValue << " at dofNoLocal: " << dofNoLocal << ", component: " << aComponent;
        }
      }  // aComponent
    }  // aDof


    // add entries in result vector for pressure
    // loop over indices of unknows (aDof,aComponent)
    for (int aDof = 0; aDof < nPressureDofsPerElement; aDof++)      // L
    {
      // get integrated value
      double integratedValue = integratedValuesPressure[aDof];

      // get local dof no, aDof is the dof within the element, dofNoLocal is the dof within the local subdomain
      dof_no_t dofNoLocal = pressureDofNosLocal[aDof];

      // set value of result vector
      if (useNestedMat_)
      {
        ierr = VecSetValue(subvectorsResidual_[3], dofNoLocal, integratedValue, ADD_VALUES); CHKERRV(ierr);
      }
      else
      {
        combinedVecResidual_->setValue(3, dofNoLocal, integratedValue, ADD_VALUES);
        VLOG(1) << "p: set value " << integratedValue << " at dofNoLocal: " << dofNoLocal;
      }
    }
  }  // elementNoLocal

  // assemble result vector
  if (useNestedMat_)
  {
    for (int componentNo = 0; componentNo < 4; componentNo++)    // a
    {
      ierr = VecAssemblyBegin(subvectorsResidual_[componentNo]); CHKERRV(ierr);
    }
    for (int componentNo = 0; componentNo < 4; componentNo++)    // a
    {
      ierr = VecAssemblyEnd(subvectorsResidual_[componentNo]); CHKERRV(ierr);
    }
  }
  else
  {
    combinedVecResidual_->finishGhostManipulation();
  }

  // if useNestedMat_: here, subvectorsResidual_[0-2] contains δW_int, subvectorsResidual_[3] contains int_Ω (J-1)*Ψ dV
  // else solverVariableResidual_ contains δW_int

  // output values for debugging
  if (outputValues)
  {
    //LOG(DEBUG) << "input u: " << *this->data_.displacements();
    //LOG(DEBUG) << "input p: " << *this->data_.pressure();

    if (useNestedMat_)
    {
      for (int componentNo = 0; componentNo < 3; componentNo++)    // a
      {
        LOG(DEBUG) << "δW_int_" << std::string(1, char('x'+componentNo)) << ": " << PetscUtility::getStringVector(subvectorsResidual_[componentNo]);
      }
      for (int componentNo = 0; componentNo < 3; componentNo++)    // a
      {
        LOG(DEBUG) << "δW_ext_" << std::string(1, char('x'+componentNo)) << ": " << PetscUtility::getStringVector(neumannBoundaryConditions_->rhs()->valuesGlobal(componentNo));
      }
    }
    else
    {
      LOG(DEBUG) << "δW_int: " << getString(solverVariableResidual_);
      LOG(DEBUG) << "δW_ext: " << getString(externalVirtualWork_);
    }
  }

  // substract external virtual work from subvectorsResidual_[0-2], only traction term,
  // δW_ext = int_∂Ω T_a phi_L dS
  if (useNestedMat_)
  {
    for (int componentNo = 0; componentNo < 3; componentNo++)    // a
    {
      ierr = VecAXPY(subvectorsResidual_[componentNo], -1, neumannBoundaryConditions_->rhs()->valuesGlobal(componentNo)); CHKERRV(ierr);
    }

    // assemble result vector
    for (int componentNo = 0; componentNo < 4; componentNo++)    // a
    {
      ierr = VecAssemblyBegin(subvectorsResidual_[componentNo]); CHKERRV(ierr);
    }
    for (int componentNo = 0; componentNo < 4; componentNo++)    // a
    {
      ierr = VecAssemblyEnd(subvectorsResidual_[componentNo]); CHKERRV(ierr);
    }
  }
  else
  {
    ierr = VecAXPY(solverVariableResidual_, -1, externalVirtualWork_); CHKERRV(ierr);
  }

  if(outputValues)
    LOG(DEBUG) << "total:   " << getString(solverVariableResidual_);

  if (useNestedMat_)
  {
    ierr = VecAssemblyBegin(solverVariableResidual_); CHKERRV(ierr);
    ierr = VecAssemblyEnd(solverVariableResidual_); CHKERRV(ierr);
  }

  // dump output vector
  // dumpVector(std::string filename, std::string format, Vec &vector, MPI_Comm mpiCommunicator, int componentNo=0, int nComponents=1);
  filename.str("");
  filename << "out/F" << std::setw(3) << std::setfill('0') << evaluationNo;
  //PetscUtility::dumpVector(filename.str(), "matlab", solverVariableResidual_, displacementsFunctionSpace->meshPartition()->mpiCommunicator());
  combinedVecResidual_->dumpGlobalNatural(filename.str());

  assert(combinedVecResidual_->currentRepresentation() == Partition::values_representation_t::representationCombinedGlobal);
  assert(combinedVecSolution_->currentRepresentation() == Partition::values_representation_t::representationCombinedGlobal);

  evaluationNo++;
}

void QuasiStaticHyperelasticitySolver::
materialComputeJacobian()
{
  // compute jacobian in submatrices_

}

Tensor2<3> QuasiStaticHyperelasticitySolver::
computeDeformationGradient(const std::array<Vec3,DisplacementsFunctionSpace::nDofsPerElement()> &displacements,
                           const Tensor2<3> &inverseJacobianMaterial,
                           const std::array<double, 3> xi
                          )
{
  // compute the deformation gradient x_i,j = δ_ij + u_i,j
  // where j is dimensionColumn and i is component of the used Vec3's


  VLOG(3) << "compute deformation gradient Fij, displacements: " << displacements;

  const int D = 3;
  Tensor2<D> deformationGradient;

  const int nDofsPerElement = DisplacementsFunctionSpace::nDofsPerElement();

  // loop over dimension, i.e. columns of deformation gradient, j
  for (int dimensionColumn = 0; dimensionColumn < 3; dimensionColumn++)
  {
    // compute du_i/dX_dimensionColumn for all i at once
    Vec3 du_dX({0});    // vector contains entries for all i, du_dXj with j = dimensionColumn

    VLOG(3) << " j = " << dimensionColumn;

    // loop over dof no., M
    for (int dofIndex = 0; dofIndex < nDofsPerElement; dofIndex++)
    {
      VLOG(3) << "  M = " << dofIndex;

      // compute dphi_dofIndex/dX_dimensionColumn
      double dphi_dX = 0;
      for (int l = 0; l < 3; l++)
      {
        VLOG(3) << "   l = " << l;
        double dphi_dxil = DisplacementsFunctionSpace::dphi_dxi(dofIndex, l, xi);
        double dxil_dX = inverseJacobianMaterial[dimensionColumn][l];     // inverseJacobianMaterial[j][l] = J_lj = dxi_l/dX_j

        VLOG(3) << "     dphi_dxil = " << dphi_dxil << ", dxil_dX = " << dxil_dX;

        // multiply dphi/dxi with dxi/dX to obtain dphi/dX
        dphi_dX += dphi_dxil * dxil_dX;
      }
      VLOG(3) << "   dphi_dX = " << dphi_dX;

      VLOG(3) << "   displ_M: " << displacements[dofIndex];
      du_dX += dphi_dX * displacements[dofIndex];   // vector-valued addition
    }
    VLOG(3) << " du_dXj: " << du_dX;

    deformationGradient[dimensionColumn] = du_dX;

    // add Kronecker delta to obtain x_i,j = delta_ij + u_i,j
    deformationGradient[dimensionColumn][dimensionColumn] += 1;
  }

  VLOG(3) << "deformationGradient: " << deformationGradient;
  return deformationGradient;
}

Tensor2<3> QuasiStaticHyperelasticitySolver::
computeRightCauchyGreenTensor(const Tensor2<3> &deformationGradient)
{
  // compute C = F^T*F where F is the deformationGradient and C is the right Cauchy-Green Tensor
  // the quantities are 3x3 tensors for the 3D case
  Tensor2<3> rightCauchyGreenTensor({std::array<double,3>({0})});

  // C_ji
  // loop over dimension, i.e. columns of right cauchy green tensor, i
  for (int dimensionColumn = 0; dimensionColumn < 3; dimensionColumn++)
  {
    // loop over row of tensor, j
    for (int dimensionRow = 0; dimensionRow < 3; dimensionRow++)
    {
      for (int k = 0; k < 3; k++)
      {
        // C_cr = C_rc += (F^T)_ck * F_kr = F_kc * F_kr
        rightCauchyGreenTensor[dimensionColumn][dimensionRow] +=
          deformationGradient[dimensionColumn][k] * deformationGradient[dimensionRow][k];
      }
    }
  }

  return rightCauchyGreenTensor;
}

std::array<double,3> QuasiStaticHyperelasticitySolver::
computeInvariants(const Tensor2<3> &rightCauchyGreen, const double rightCauchyGreenDeterminant)
{
  std::array<double,3> invariants;

  // I1 = tr(C)
  invariants[0] = 0.0;
  for (int i = 0; i < 3; i++)
  {
    invariants[0] += rightCauchyGreen[i][i];
  }

  // tr(C^2)
  double traceC2 = 0;
  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      traceC2 += rightCauchyGreen[i][j] * rightCauchyGreen[j][i];
    }
  }

  // I2 = 1/2 * (tr(C)^2 - tr(C^2))
  invariants[1] = 0.5 * (MathUtility::sqr(invariants[0]) - traceC2);

  // I3 = det(C)
  invariants[2] = rightCauchyGreenDeterminant;

  return invariants;
}

std::array<double,2> QuasiStaticHyperelasticitySolver::
computeReducedInvariants(const std::array<double,3> invariants, const double deformationGradientDeterminant)
{
  std::array<double,2> reducedInvariants;

  // 3D: Fbar = J^{-1/3}*F such that det(Fbar) = 1
  // 2D: Fbar = J^{-1/2}*F such that det(Fbar) = 1

  // 3D: Cbar = J^{-2/3}*C such that det(Cbar) = 1
  // 2D: Cbar = J^{-2/2}*C = J^{-1}*C = J^{-1}*F^T*F = J^{-1}*J^{1/2}*J^{1/2}*Fbar^T*Fbar = Fbar^T*Fbar such that det(Cbar) = 1

  double factor23 = -2./3;
  double factor43 = -4./3;

  // for 2D, not used here
#if 0
  factor23 = -1./2;
  factor43 = -2./2;
#endif

  reducedInvariants[0] = pow(deformationGradientDeterminant, factor23) * invariants[0];
  reducedInvariants[1] = pow(deformationGradientDeterminant, factor43) * invariants[1];

  return reducedInvariants;
}

Tensor2<3> QuasiStaticHyperelasticitySolver::
computePK2Stress(const double pressure,                             //< [in] pressure value p
                 const Tensor2<3> &rightCauchyGreen,                //< [in] C
                 const Tensor2<3> &inverseRightCauchyGreen,         //< [in] C^{-1}
                 const std::array<double,2> reducedInvariants,      //< [in] the reduced invariants Ibar_1, Ibar_2
                 const double deformationGradientDeterminant        //< [in] J = det(F)
                )
{
  // compute the PK2 stress tensor as S=2*dPsi/dC
  // for explanation see pdf document
  typedef Equation::Static::MooneyRivlinIncompressible3D Term;

  auto dPsi_dIbar1Expression = SEMT::deriv_t(Term::strainEnergyDensityFunctionIsochoric, Term::Ibar1);
  auto dPsi_dIbar2Expression = SEMT::deriv_t(Term::strainEnergyDensityFunctionIsochoric, Term::Ibar2);

  std::vector<double> reducedInvariantsVector(reducedInvariants.begin(), reducedInvariants.end());

  const double dPsi_dIbar1 = dPsi_dIbar1Expression.apply(reducedInvariantsVector);
  const double dPsi_dIbar2 = dPsi_dIbar2Expression.apply(reducedInvariantsVector);

  const double Ibar1 = reducedInvariants[0];
  const double J = deformationGradientDeterminant;

  double factor23 = -2./3;
  double factor1 = 2*(dPsi_dIbar1 + Ibar1*dPsi_dIbar2);
  double factor2 = -2*dPsi_dIbar2;
  double factorJ23 = pow(J, factor23);

  Tensor2<3> fictitiousPK2Stress;          //< [out] Sbar, the fictitious 2nd Piola-Kirchhoff stress tensor
  //Tensor2<3> pk2StressIsochoric            //< [out] S_iso, the isochoric part of the 2nd Piola-Kirchhoff stress tensor
  Tensor2<3> pK2Stress;

  // compute fictitiousPK2Stress:
  // Sbar = factor1*I + factor2*Cbar, Cbar = J^{-2/3}*C

  // row index
  for (int i=0; i<3; i++)
  {
    // column index
    for (int j=0; j<3; j++)
    {
      int delta_ij = (i == j? 1 : 0);
      const double cBar = factorJ23 * rightCauchyGreen[j][i];
      fictitiousPK2Stress[j][i] = factor1 * delta_ij + factor2 * cBar;

      //if (i == j)
      //  LOG(DEBUG) << "  C_" << i << j << " = " << rightCauchyGreen[j][i] << ", factorJ23 = " << factorJ23 << ", SBar_" << i << j << ": " << fictitiousPK2Stress[j][i]
      //    << " = " << factor1 << " * " << delta_ij << " + " << factor2 << "*" << cBar ;
    }
  }

  // Holzapfel p.234
  // S = S_vol + S_iso
  // S_vol = J*p*C^-1
  // S_iso = 2*dPsi_iso/dC = J^(-2/3) P : Sbar   (P = II - 1/3(C^{-1} dyad C): Holzapfel p.229)
  // pSbar = P : Sbar, P: 4th order tensor, Sbar: 2nd order tensor, note: A:B = A_ijkl*B_kl*e_i dyad e_j

  // compute S = S_vol + S_iso
  // row index
  for (int i=0; i<3; i++)     // alternative indices: A
  {
    // column index
    for (int j=0; j<3; j++)     // alternative indices: B
    {
      // volumetric stress
      const double sVol = J * pressure * inverseRightCauchyGreen[j][i];         // S_vol = J * p * C^{-1}_AB

      // compute P : Sbar
      double pSbar = 0;
      //double ccs = 0;  // (C^-1 dyad C)*Sbar
      //double cSbar = 0;  // C:Sbar

      // row index
      for (int k=0; k<3; k++)        // alternative indices: C
      {
        const int delta_ik = (i == k? 1 : 0);
        const int delta_jk = (j == k? 1 : 0);

        // column index
        for (int l=0; l<3; l++)            // alternative indices: D
        {
          const int delta_il = (i == l? 1 : 0);
          const int delta_jl = (j == l? 1 : 0);
          const int Ii = 0.5 * (delta_ik * delta_jl + delta_il * delta_jk);       // II = (δ_AC*δ_BD + δ_AD*δ_BC) / 2

          const double Cc = inverseRightCauchyGreen[j][i] * rightCauchyGreen[l][k];     // CC = C^{-1}_AB * C_CD
          const double Pp = (Ii - 1./3 * Cc);

          pSbar += Pp * fictitiousPK2Stress[l][k];

          //cSbar += rightCauchyGreen[l][k] * fictitiousPK2Stress[l][k];
          //ccs += Cc * fictitiousPK2Stress[l][k];
        }
      }

      // isochoric stress
      const double sIso = factorJ23 * pSbar;
      //pk2StressIsochoric[j][i] = sIso;

      // total stress is sum of volumetric and isochoric part, sVol = J*p*C^{-1}_AB, sIso = j^{-2/3}*(Ii-1/3*Cc)*Sbar
      pK2Stress[j][i] = sVol + sIso;
      //VLOG(2) << "set pk2Stress_" << i << j << " = " << pK2Stress[j][i];

      //if (i == j)
      //  LOG(DEBUG) << "  ccs: " << ccs << " C:Sbar: " << cSbar << ", factorJ23: " << factorJ23 << ", Svol_" << i << j << " = " << sVol << ", Siso_" << i << j << " = " << sIso << ", S = " << pK2Stress[j][i];
    }
  }

  //for debugging, check symmetry of PK2 stress and if it is correct
  const double errorTolerance = 1e-14;
  if (VLOG_IS_ON(2))
  {
    bool pK2IsSymmetric = true;
    for (int a=0; a<3; a++)
    {
      for (int b=0; b<3; b++)
      {
        if (fabs(pK2Stress[a][b] - pK2Stress[b][a]) > errorTolerance)
        {
          LOG(ERROR) << "pK2Stress[" <<a<< "][" <<b<< "] != pK2Stress[" <<b<< "][" <<a<< "] (" <<pK2Stress[b][a]<< " != " <<pK2Stress[a][b]<< ")";
          pK2IsSymmetric = false;
        }
      }
    }
    if (pK2IsSymmetric)
      VLOG(2) << "PK2 stress tensor is symmetric!";

    // check if PK2 is the same as from explicit formula for Mooney-Rivlin (p.249)

    // explicit formula in Holzapfel p.249

    double factor23 = -2./3;
    const double factorJ23 = pow(deformationGradientDeterminant, factor23);

    const double c0 = SEMT::Parameter<0>::get_value();
    const double c1 = SEMT::Parameter<1>::get_value();
    //const double c0 = PARAM(0).get_value();   //< material parameter
    //const double c1 = PARAM(1).get_value();   //< material parameter

    const double Ibar1 = reducedInvariants[0];

    const double gamma1 = 2*(c0 + c1*Ibar1);
    const double gamma2 = -2*c1;

    const double errorTolerance = 1e-14;
    bool mismatch = false;
    for (int a = 0; a < 3; a++)
    {
      for (int b = 0; b < 3; b++)
      {
        const int delta_ab = (a == b? 1 : 0);
        double cBar = factorJ23 * rightCauchyGreen[b][a];
        double sBar = gamma1*delta_ab + gamma2*cBar;

        if (fabs(sBar - fictitiousPK2Stress[b][a]) > errorTolerance)
        {
          LOG(ERROR) << "mismatch in Sbar_" << a << b << ": derived: " << fictitiousPK2Stress << ", explicit formula: " << sBar;
          mismatch = true;
        }
      }
    }
    if (!mismatch)
      VLOG(2) << "Sbar is correct!";
  }

  return pK2Stress;
}

void QuasiStaticHyperelasticitySolver::
computePK2StressField()
{
  LOG(TRACE) << "computePK2StressField";

  // get pointer to function space
  std::shared_ptr<DisplacementsFunctionSpace> displacementsFunctionSpace = this->data_.displacementsFunctionSpace();
  std::shared_ptr<PressureFunctionSpace> pressureFunctionSpace = this->data_.pressureFunctionSpace();

  const int D = 3;  // dimension
  const int nDisplacementsDofsPerElement = DisplacementsFunctionSpace::nDofsPerElement();
  const int nPressureDofsPerElement = PressureFunctionSpace::nDofsPerElement();
  const int nElementsLocal = displacementsFunctionSpace->nElementsLocal();

  // loop over elements
  for (int elementNoLocal = 0; elementNoLocal < nElementsLocal; elementNoLocal++)
  {
    LOG(DEBUG) << "el " << elementNoLocal;
    // get geometry field of reference configuration
    std::array<Vec3,nDisplacementsDofsPerElement> geometryReferenceValues;
    this->data_.geometryReference()->getElementValues(elementNoLocal, geometryReferenceValues);

    // get displacements field values for element
    std::array<Vec3,nDisplacementsDofsPerElement> displacementsValues;
    this->data_.displacements()->getElementValues(elementNoLocal, displacementsValues);

    std::array<double,nPressureDofsPerElement> pressureValuesCurrentElement;
    this->data_.pressure()->getElementValues(elementNoLocal, pressureValuesCurrentElement);

    // get indices of element-local dofs
    std::array<dof_no_t,27> dofNosLocal = this->displacementsFunctionSpace_->getElementDofNosLocal(elementNoLocal);

    // loop over nodes of this element
    for (int elementalNodeNo = 0; elementalNodeNo < 27; elementalNodeNo++)
    {
      dof_no_t dofNoLocal = dofNosLocal[elementalNodeNo];

      // get parameter values of current sampling point
      int indexX = elementalNodeNo % 3;
      int indexY = int((elementalNodeNo % 9) / 3);
      int indexZ = int(elementalNodeNo / 9);
      Vec3 xi({});

      switch (indexX)
      {
      case 0:
        xi[0] = 0.0;
        break;
      case 1:
        xi[0] = 0.5;
        break;
      case 2:
        xi[0] = 1.0;
        break;
      };
      switch (indexY)
      {
      case 0:
        xi[1] = 0.0;
        break;
      case 1:
        xi[1] = 0.5;
        break;
      case 2:
        xi[1] = 1.0;
        break;
      };
      switch (indexZ)
      {
      case 0:
        xi[2] = 0.0;
        break;
      case 1:
        xi[2] = 0.5;
        break;
      case 2:
        xi[2] = 1.0;
        break;
      };

      // compute the 3x3 jacobian of the parameter space to world space mapping
      Tensor2<D> jacobianMaterial = DisplacementsFunctionSpace::computeJacobian(geometryReferenceValues, xi);
      double jacobianDeterminant;
      Tensor2<D> inverseJacobianMaterial = MathUtility::computeInverse<D>(jacobianMaterial, jacobianDeterminant);

      // jacobianMaterial[columnIdx][rowIdx] = dX_rowIdx/dxi_columnIdx
      // inverseJacobianMaterial[columnIdx][rowIdx] = dxi_rowIdx/dX_columnIdx because of inverse function theorem

      // F
      Tensor2<D> deformationGradient = this->computeDeformationGradient(displacementsValues, inverseJacobianMaterial, xi);
      double deformationGradientDeterminant = MathUtility::computeDeterminant<D>(deformationGradient);  // J

      Tensor2<D> rightCauchyGreen = this->computeRightCauchyGreenTensor(deformationGradient);  // C = F^T*F

      double rightCauchyGreenDeterminant;   // J^2
      Tensor2<D> inverseRightCauchyGreen = MathUtility::computeSymmetricInverse<D>(rightCauchyGreen, rightCauchyGreenDeterminant);  // C^-1

      // invariants
      std::array<double,3> invariants = this->computeInvariants(rightCauchyGreen, rightCauchyGreenDeterminant);  // I_1, I_2, I_3
      std::array<double,2> reducedInvariants = this->computeReducedInvariants(invariants, deformationGradientDeterminant); // Ibar_1, Ibar_2

      // pressure is the separately interpolated pressure for mixed formulation
      double pressure = pressureFunctionSpace->interpolateValueInElement(pressureValuesCurrentElement, xi);

      // Pk2 stress tensor S = S_vol + S_iso (p.234)
      //! compute 2nd Piola-Kirchhoff stress tensor S = 2*dPsi/dC and the fictitious PK2 Stress Sbar
      Tensor2<D> pK2Stress = this->computePK2Stress(pressure, rightCauchyGreen, inverseRightCauchyGreen, reducedInvariants, deformationGradientDeterminant);

      std::array<double,6> valuesInVoigtNotation({pK2Stress[0][0], pK2Stress[1][1], pK2Stress[2][2], pK2Stress[0][1], pK2Stress[1][2], pK2Stress[0][2]});

      LOG(DEBUG) << "node " << dofNoLocal << " pk2: " << valuesInVoigtNotation;
      this->data_.pK2Stress()->setValue(dofNoLocal, valuesInVoigtNotation, INSERT_VALUES);
    }
  }
  this->data_.pK2Stress()->finishGhostManipulation();
}

} // namespace
