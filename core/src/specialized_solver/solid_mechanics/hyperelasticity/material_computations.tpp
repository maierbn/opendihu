#include "specialized_solver/solid_mechanics/hyperelasticity/hyperelasticity_solver.h"

#include <Python.h>  // has to be the first included header
#include <array>
#include <Vc/Vc>

#include "equation/mooney_rivlin_incompressible.h"

namespace SpatialDiscretization
{


template<typename Term,typename MeshType, int nDisplacementComponents>
bool HyperelasticitySolver<Term,MeshType,nDisplacementComponents>::
materialComputeInternalVirtualWork(bool communicateGhosts)
{
  // compute Wint in solverVariableResidual_
  //  output is solverVariableResidual_, a normal Vec, no Dirichlet BC dofs, also accessible by combinedVecResidual_
  //  input is solverVariableSolution_, a normal Vec, the same values have already been assigned to this->data_.displacements() and this->data_.pressure()

  // -----------------------------
  // compute internal virtual work
  // δW_int = int_Ω 1/2 S_AB (F_aB phi_L,A + F_aA phi_L,B) dV

  const bool outputFiles = false;

  if (communicateGhosts)
  {
    // assert that data representation is global
    if (combinedVecResidual_->currentRepresentation() != Partition::values_representation_t::representationCombinedGlobal)
      LOG(ERROR) << "Representation is " << combinedVecResidual_->getCurrentRepresentationString();
    assert(combinedVecResidual_->currentRepresentation() == Partition::values_representation_t::representationCombinedGlobal);
    assert(combinedVecSolution_->currentRepresentation() == Partition::values_representation_t::representationCombinedGlobal);
  }
  VLOG(1) << "materialComputeInternalVirtualWork";

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
  typedef std::array<double_v_t, nUnknowsPerElement> EvaluationsDisplacementsType;
  std::array<EvaluationsDisplacementsType, QuadratureDD::numberEvaluations()> evaluationsArrayDisplacements{};

  typedef std::array<double_v_t, nPressureDofsPerElement> EvaluationsPressureType;
  std::array<EvaluationsPressureType, QuadratureDD::numberEvaluations()> evaluationsArrayPressure{};

  // setup arrays used for integration
  std::array<Vec3, QuadratureDD::numberEvaluations()> samplingPoints = QuadratureDD::samplingPoints();

  // set values to zero
  if (communicateGhosts)
  {
    combinedVecResidual_->zeroEntries();
    combinedVecResidual_->startGhostManipulation();
  }

  static int evaluationNo = 0;  // counter how often this function was called

  if (outputFiles)
  {
    // dump input vector
    // dumpVector(std::string filename, std::string format, Vec &vector, MPI_Comm mpiCommunicator, int componentNo=0, int nComponents=1);
    std::stringstream filename;
    filename << "out/x" << std::setw(3) << std::setfill('0') << evaluationNo++;
    //PetscUtility::dumpVector(filename.str(), "matlab", solverVariableSolution_, displacementsFunctionSpace->meshPartition()->mpiCommunicator());
    combinedVecSolution_->dumpGlobalNatural(filename.str());
  }

  // loop over elements, always 4 elements at once using the vectorized functions
  for (int elementNoLocal = 0; elementNoLocal < nElementsLocal; elementNoLocal += nVcComponents)
  {

#ifdef USE_VECTORIZED_FE_MATRIX_ASSEMBLY
    // get indices of elementNos that should be handled in the current iterations,
    // this is, e.g.
    //    [10,11,12,13,-1,-1,-1,-1] (if nVcComponents==4 and nElementsLocal > 13)
    // or [10,11,12,-1,-1,-1,-1,-1] (if nVcComponents==4 and nElementsLocal == 13)


    dof_no_v_t elementNoLocalv([elementNoLocal, nElementsLocal](int i)
    {
      return (i >= nVcComponents || elementNoLocal+i >= nElementsLocal? -1: elementNoLocal+i);
    });

    // here, elementNoLocalv is the list of indices of the current iteration, e.g. [10,11,12,13,-1,-1,-1,-1]
    // elementNoLocal is the first entry of elementNoLocalv
#else
    int elementNoLocalv = elementNoLocal;
#endif

    // get geometry field of reference configuration
    std::array<Vec3_v_t,nDisplacementsDofsPerElement> geometryReferenceValues;
    this->data_.geometryReference()->getElementValues(elementNoLocalv, geometryReferenceValues);

    // get displacements field values for element
    std::array<Vec3_v_t,nDisplacementsDofsPerElement> displacementsValues;
    this->data_.displacements()->getElementValues(elementNoLocalv, displacementsValues);

    if (VLOG_IS_ON(1))
    {
      global_no_t elementNoGlobal = displacementsFunctionSpace->meshPartition()->getElementNoGlobalNatural(elementNoLocal);
      VLOG(1) << "elementNoLocal " << elementNoLocal << ", displacementsValues: " << displacementsValues;
      VLOG(1) << "elementNoGlobal " << elementNoGlobal << ", geometryReferenceValues: " << geometryReferenceValues;
    }

    std::array<double_v_t,nPressureDofsPerElement> pressureValuesCurrentElement;
    this->data_.pressure()->getElementValues(elementNoLocalv, pressureValuesCurrentElement);

    std::array<Vec3_v_t,nDisplacementsDofsPerElement> elementalDirectionValues;
    this->data_.fiberDirection()->getElementValues(elementNoLocalv, elementalDirectionValues);

    std::array<VecD_v_t<6>,nDisplacementsDofsPerElement> activePK2StressValues;
    if (Term::usesActiveStress)
    {
      this->data_.activePK2Stress()->getElementValues(elementNoLocalv, activePK2StressValues);
    }

    // loop over integration points (e.g. gauss points) for displacements field
    for (unsigned int samplingPointIndex = 0; samplingPointIndex < samplingPoints.size(); samplingPointIndex++)
    {
      // get parameter values of current sampling point
      Vec3 xi = samplingPoints[samplingPointIndex];

      // compute the 3x3 jacobian of the parameter space to world space mapping
      Tensor2_v_t<D> jacobianMaterial = DisplacementsFunctionSpace::computeJacobian(geometryReferenceValues, xi);
      double_v_t jacobianDeterminant;
      Tensor2_v_t<D> inverseJacobianMaterial = MathUtility::computeInverse(jacobianMaterial, jacobianDeterminant);

      // jacobianMaterial[columnIdx][rowIdx] = dX_rowIdx/dxi_columnIdx
      // inverseJacobianMaterial[columnIdx][rowIdx] = dxi_rowIdx/dX_columnIdx because of inverse function theorem

      // get the factor in the integral that arises from the change in integration domain from world to parameter space
      double_v_t integrationFactor = MathUtility::abs(jacobianDeterminant); // MathUtility::computeIntegrationFactor(jacobianMaterial);

      // F
      Tensor2_v_t<D> deformationGradient = this->computeDeformationGradient(displacementsValues, inverseJacobianMaterial, xi);
      double_v_t deformationGradientDeterminant = MathUtility::computeDeterminant(deformationGradient);  // J

      Tensor2_v_t<D> rightCauchyGreen = this->computeRightCauchyGreenTensor(deformationGradient);  // C = F^T*F

      double_v_t rightCauchyGreenDeterminant;   // J^2
      Tensor2_v_t<D> inverseRightCauchyGreen = MathUtility::computeSymmetricInverse(rightCauchyGreen, rightCauchyGreenDeterminant);  // C^-1

      // fiber direction
      Vec3_v_t fiberDirection = displacementsFunctionSpace->template interpolateValueInElement<3>(elementalDirectionValues, xi);

#ifndef NDEBUG
      if (Term::usesFiberDirection)
      {
        if (fabs(MathUtility::norm<3>(fiberDirection) - 1) > 1e-3)
          LOG(FATAL) << "fiberDirecton " << fiberDirection << " is not normalized, elementalDirectionValues:" << elementalDirectionValues;
      }
#endif

      // invariants
      std::array<double_v_t,5> invariants = this->computeInvariants(rightCauchyGreen, rightCauchyGreenDeterminant, fiberDirection);  // I_1, I_2, I_3, I_4, I_5
      std::array<double_v_t,5> reducedInvariants = this->computeReducedInvariants(invariants, deformationGradientDeterminant); // Ibar_1, Ibar_2, Ibar_4, Ibar_5

      // pressure is the separately interpolated pressure for mixed formulation
      double_v_t pressure = pressureFunctionSpace->interpolateValueInElement(pressureValuesCurrentElement, xi);

      // Pk2 stress tensor S = S_vol + S_iso (p.234)
      //! compute 2nd Piola-Kirchhoff stress tensor S = 2*dPsi/dC and the fictitious PK2 Stress Sbar
      Tensor2_v_t<D> fictitiousPK2Stress;   // Sbar
      Tensor2_v_t<D> pk2StressIsochoric;    // S_iso
      Tensor2_v_t<D> pK2Stress = this->computePK2Stress(pressure, rightCauchyGreen, inverseRightCauchyGreen, invariants, reducedInvariants,
                                                        deformationGradientDeterminant, fiberDirection,
                                                        fictitiousPK2Stress, pk2StressIsochoric);

      // add active stress contribution if this material has this
      if (Term::usesActiveStress)
      {
        VecD_v_t<6> activePK2StressInVoigtNotation = displacementsFunctionSpace->template interpolateValueInElement<6>(activePK2StressValues, xi);

        pK2Stress[0][0] += activePK2StressInVoigtNotation[0];
        pK2Stress[1][1] += activePK2StressInVoigtNotation[1];
        pK2Stress[2][2] += activePK2StressInVoigtNotation[2];
        pK2Stress[0][1] += activePK2StressInVoigtNotation[3];
        pK2Stress[1][0] += activePK2StressInVoigtNotation[3];
        pK2Stress[1][2] += activePK2StressInVoigtNotation[4];
        pK2Stress[2][1] += activePK2StressInVoigtNotation[4];
        pK2Stress[0][2] += activePK2StressInVoigtNotation[5];
        pK2Stress[2][0] += activePK2StressInVoigtNotation[5];
      }

      // call debugging methods, currently disabled
      this->materialTesting(pressure, rightCauchyGreen, inverseRightCauchyGreen, reducedInvariants, deformationGradientDeterminant, fiberDirection,
        fictitiousPK2Stress, pk2StressIsochoric
      );

      std::array<Vec3,nDisplacementsDofsPerElement> gradPhi = displacementsFunctionSpace->getGradPhi(xi);
      // (column-major storage) gradPhi[L][a] = dphi_L / dxi_a
      // gradPhi[column][row] = gradPhi[dofIndex][i] = dphi_dofIndex/dxi_i, columnIdx = dofIndex, rowIdx = which direction

      if (VLOG_IS_ON(2))
      {
        global_no_t elementNoGlobal = displacementsFunctionSpace->meshPartition()->getElementNoGlobalNatural(elementNoLocal);

        VLOG(2) << "";
        VLOG(2) << "element local " << elementNoLocal << " global " << elementNoGlobal << " xi: " << xi;
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
        VLOG(2) << "  PK2Stress: S=" << pK2Stress;
        VLOG(2) << "  gradPhi: " << gradPhi;
      }

      VLOG(1) << "  sampling point " << samplingPointIndex << "/" << samplingPoints.size() << ", xi: " << xi << ", J: " << deformationGradientDeterminant << ", p: " << pressure << ", S11: " << pK2Stress[0][0];

      if (samplingPointIndex == 0 && D == 3)
        VLOG(1) << " F11: " << deformationGradient[0][0] << ", F22,F33: " << deformationGradient[1][1] << "," << deformationGradient[2][2]
          << ", F12,F13,F23: " << deformationGradient[1][0] << "," << deformationGradient[2][0] << "," << deformationGradient[2][1]
          << ", J: " << deformationGradientDeterminant << ", p: " << pressure;

      /*if (VLOG_IS_ON(2))
      {
        Tensor2<D> greenLangrangeStrain = this->computeGreenLagrangeStrain(rightCauchyGreen);
        VLOG(2) << "  strain E=" << greenLangrangeStrain;
      }*/

      if (Vc::any_of(deformationGradientDeterminant < 1e-12))   // if any entry of the deformation gradient is negative
      {
        LOG(WARNING) << "Deformation gradient " << deformationGradient << " has zero or negative determinant " << deformationGradientDeterminant
          << std::endl << "Geometry values in element " << elementNoLocal << ": " << geometryReferenceValues << std::endl
          << "Displacements at xi " << xi << ": " << displacementsValues;
        lastSolveSucceeded_ = false;
      }

      // loop over basis functions and evaluate integrand at xi for displacement part (δW_int - δW_ext)
      for (int aDof = 0; aDof < nDisplacementsDofsPerElement; aDof++)           // index over dofs, each dof has D components, L in derivation
      {
        for (int aComponent = 0; aComponent < D; aComponent++)     // lower-case a in derivation, index over displacements components
        {
          // compute index of degree of freedom and component (result vector index)
          const int i = D*aDof + aComponent;

          // compute result[i]
          double_v_t integrand = 0.0;
          for (int aInternal = 0; aInternal < D; aInternal++)     // capital A in derivation
          {
            for (int bInternal = 0; bInternal < D; bInternal++)     // capital B in derivation
            {
              const double_v_t faB = deformationGradient[bInternal][aComponent];
              const double_v_t faA = deformationGradient[aInternal][aComponent];

              // ----------------------------
              // compute derivatives of phi
              // note that dphi^L_a = dphi^L, i.e. dphi^L_{b,A} = dphi^L_{c,A} = dphi^L_{,A}
              double_v_t dphiL_dXA = 0.0;
              double_v_t dphiL_dXB = 0.0;

              // helper index k for multiplication with inverse Jacobian
              for (int k = 0; k < D; k++)
              {
                // compute dphiL/dXA from dphiL/dxik and dxik/dXA
                const double dphiL_dxik = gradPhi[aDof][k];    // dphi_L/dxik
                const double_v_t dxik_dXA = inverseJacobianMaterial[aInternal][k];  // inverseJacobianMaterial[A][k] = J^{-1}_kA = dxi_k/dX_A

                dphiL_dXA += dphiL_dxik * dxik_dXA;

                // compute dphiL/dXB from dphiL/dxik and dxik/dXB
                const double_v_t dxik_dXB = inverseJacobianMaterial[bInternal][k];  // inverseJacobianMaterial[B][k] = J^{-1}_kB = dxi_k/dX_B

                dphiL_dXB += dphiL_dxik * dxik_dXB;
              }

              integrand += 1./2. * pK2Stress[bInternal][aInternal] * (faB*dphiL_dXA + faA*dphiL_dXB);

            }  // B, bInternal
          }  // A, aInternal

          VLOG(2) << "   (L,a)=(" << aDof << "," << aComponent << "), integrand: " << integrand;

          // store integrand in evaluations array
          evaluationsArrayDisplacements[samplingPointIndex][i] = integrand * integrationFactor;

        }  // a
      }  // L

      // for the incompressible formulation, also integrate the incompressibility constraint
      if (Term::isIncompressible)
      {
        // loop over basis functions and evaluate integrand at xi for pressure part ((J-1)*psi)
        for (int dofIndex = 0; dofIndex < nPressureDofsPerElement; dofIndex++)           // index over dofs in element, L in derivation
        {
          const double phiL = pressureFunctionSpace->phi(dofIndex, xi);
          const double_v_t integrand = (deformationGradientDeterminant - 1.0) * phiL;     // (J-1) * phi_L

          // store integrand in evaluations array
          evaluationsArrayPressure[samplingPointIndex][dofIndex] = integrand * integrationFactor;
        }  // L
      }

    }  // function evaluations

    // integrate all values for result vector entries at once
    EvaluationsDisplacementsType integratedValuesDisplacements = QuadratureDD::computeIntegral(evaluationsArrayDisplacements);

    EvaluationsPressureType integratedValuesPressure;
    if (Term::isIncompressible)
    {
      integratedValuesPressure = QuadratureDD::computeIntegral(evaluationsArrayPressure);
    }

    // get indices of element-local dofs
    std::array<dof_no_v_t,nDisplacementsDofsPerElement> dofNosLocal = displacementsFunctionSpace->getElementDofNosLocal(elementNoLocalv);
    std::array<dof_no_v_t,nPressureDofsPerElement> pressureDofNosLocal = pressureFunctionSpace->getElementDofNosLocal(elementNoLocalv);

    VLOG(2) << "  element " << elementNoLocalv << " has dofs " << dofNosLocal;

    // add entries in result vector for displacements
    // loop over indices of unknows (aDof,aComponent)
    for (int aDof = 0; aDof < nDisplacementsDofsPerElement; aDof++)      // L, this is the dof within the element
    {
      for (int aComponent = 0; aComponent < D; aComponent++)    // a
      {
        // compute index of degree of freedom and component (matrix row index)
        const int i = D*aDof + aComponent;

        // integrate value and set entry
        double_v_t integratedValue = integratedValuesDisplacements[i];

        // get local dof no, aDof is the dof within the element, dofNoLocal is the dof within the local subdomain
        dof_no_v_t dofNoLocal = dofNosLocal[aDof];

        VLOG(1) << "  result vector (L,a)=(" <<aDof<< "," <<aComponent<< "), " <<i<< ", dof " << dofNosLocal[aDof] << " all elemental dofs: "  << dofNosLocal
          << ", integrated value: " <<integratedValue;

        combinedVecResidual_->setValue(aComponent, dofNoLocal, integratedValue, ADD_VALUES);
        VLOG(1) << "u: set value " << integratedValue << " at dofNoLocal: " << dofNoLocal << ", component: " << aComponent;
      }  // aComponent
    }  // aDof

    // only for the incompressible formulation, also integrate the incompressibility constraint
    // add entries in result vector for pressure
    if (Term::isIncompressible)
    {
      // loop over indices of unknows (aDof,aComponent)
      for (int aDof = 0; aDof < nPressureDofsPerElement; aDof++)      // L
      {
        // get integrated value
        double_v_t integratedValue = integratedValuesPressure[aDof];

        // get local dof no, aDof is the dof within the element, dofNoLocal is the dof within the local subdomain
        dof_no_v_t dofNoLocal = pressureDofNosLocal[aDof];

        // set value of result vector
        const int pressureDofNo = nDisplacementComponents;  // 3 or 6, depending if static or dynamic problem
        combinedVecResidual_->setValue(pressureDofNo, dofNoLocal, integratedValue, ADD_VALUES);

        //if (VLOG_IS_ON(1))
        //{
        //  global_no_t dofNoGlobalPetsc = pressureFunctionSpace->meshPartition()->getDofNoGlobalPetsc(dofNoLocal);
        //  VLOG(1) << "p: add value " << integratedValue << " at dofNoGlobalPetsc: " << dofNoGlobalPetsc;
        //}
      }
    }  // elementNoLocal
  }

  // assemble result vector
  if (communicateGhosts)
  {
    combinedVecResidual_->finishGhostManipulation();
  }
  //combinedVecResidual_->startGhostManipulation();
  //combinedVecResidual_->zeroGhostBuffer();
  //combinedVecResidual_->finishGhostManipulation();

  if (!lastSolveSucceeded_)
  {
    // return false means the computation was not successful
    return false;
  }

  // now, solverVariableResidual_, which is the globalValues() of combinedVecResidual_, contains δW_int
  // computation was successful
  return true;
}

template<typename Term,typename MeshType, int nDisplacementComponents>
bool HyperelasticitySolver<Term,MeshType,nDisplacementComponents>::
materialComputeResidual(double loadFactor)
{
  // This computes the residual, i.e. the nonlinear function to be solved.
  // Compute Wint - Wext in variable solverVariableResidual_.
  //  output is solverVariableResidual_, a normal Vec, no Dirichlet BC dofs, also accessible by combinedVecResidual_
  //  input is solverVariableSolution_, a normal Vec, the same values have already been assigned to this->data_.displacements() and this->data_.pressure() (!)
  // before this method, values of u, v and p get stored to the data object by setUVP(solverVariableSolution_);

  LOG(DEBUG) << "materialComputeResidual(loadFactor=" << loadFactor << ")";

  const bool outputValues = false;
  const bool outputFiles = false;
  if (outputValues)
  {
    LOG(DEBUG) << "input: " << getString(solverVariableSolution_);
  }

  // dump input vector to file
  if (outputFiles)
  {
    static int evaluationNo = 0;  // counter how often this function was called

    // dumpVector(std::string filename, std::string format, Vec &vector, MPI_Comm mpiCommunicator, int componentNo=0, int nComponents=1);
    std::stringstream filename;
    filename << "out/x" << std::setw(3) << std::setfill('0') << evaluationNo++ << "r" << DihuContext::nRanksCommWorld();

    combinedVecSolution_->dumpGlobalNatural(filename.str());
  }

  // prepare combinedVecResidual_ where internal virtual work will be computed
  combinedVecResidual_->zeroEntries();
  combinedVecResidual_->startGhostManipulation();

  bool successful = materialComputeInternalVirtualWork(false);    // compute without communicating ghost values, because startGhostManipulation has been called

  // if there was a negative jacobian, exit this computation
  if (!successful)
  {
    combinedVecResidual_->finishGhostManipulation();     // communicate and add up values in ghost buffers
    return false;
  }

  // now, solverVariableResidual_, which is the globalValues() of combinedVecResidual_, contains δW_int
  // also the pressure equation residual has been set at the last component

  // output values for debugging
  if (outputValues)
  {
    //LOG(DEBUG) << "input u: " << *this->data_.displacements();
    //LOG(DEBUG) << "input p: " << *this->data_.pressure();

    LOG(DEBUG) << "δW_int: " << getString(solverVariableResidual_);
    LOG(DEBUG) << "δW_ext: " << getString(externalVirtualWorkDead_);
  }

  // for static case: F = δW_int - δW_ext

  if (nDisplacementComponents == 3)
  {
    combinedVecResidual_->finishGhostManipulation();     // communicate and add up values in ghost buffers

    // compute F = δW_int - δW_ext,
    // δW_ext = int_∂Ω T_a phi_L dS was precomputed in initialize (materialComputeExternalVirtualWorkDead()), in variable externalVirtualWorkDead_
    // for static case, externalVirtualWorkDead_ = externalVirtualWorkDead_
    PetscErrorCode ierr;
    ierr = VecAXPY(solverVariableResidual_, -loadFactor, externalVirtualWorkDead_);
    if (ierr)
      return false;

    if(outputValues)
      LOG(DEBUG) << "static problem, total F = δW_int - δW_ext:" << getString(solverVariableResidual_);

  }
  else if (nDisplacementComponents == 6)
  {
    // for dynamic case, add acceleration term to residual

    // add acceleration term (int_Ω rho_0 (v^(n+1) - v^(n)) / dt * phi^L * phi^M * δu^M dV) to solverVariableResidual_
    // also add the velocity equation in the velocity slot
    materialAddAccelerationTermAndVelocityEquation(false);

    // the incompressibility equation has been added in the pressure slot by materialComputeInternalVirtualWork

    combinedVecResidual_->finishGhostManipulation();     // communicate and add up values in ghost buffers

    // compute F = δW_int - δW_ext,dead + accelerationTerm
    // δW_ext = int_∂Ω T_a phi_L dS was precomputed in initialize, in variable externalVirtualWorkDead_
    PetscErrorCode ierr;
    ierr = VecAXPY(solverVariableResidual_, -loadFactor, externalVirtualWorkDead_);
    if (ierr)
      return false;

    if (outputFiles)
    {
      static int evaluationNo = 0;  // counter how often this function was called

      // dumpVector(std::string filename, std::string format, Vec &vector, MPI_Comm mpiCommunicator, int componentNo=0, int nComponents=1);
      std::stringstream filename;
      filename << "out/id" << std::setw(3) << std::setfill('0') << evaluationNo++ << "r" << DihuContext::nRanksCommWorld();

      combinedVecResidual_->dumpGlobalNatural(filename.str());
    }
  }

  if (outputValues)
  {
    LOG(DEBUG) << "residual: " << getString(solverVariableResidual_);
    PetscErrorCode ierr;
    PetscReal norm;
    ierr = VecNorm(solverVariableResidual_, NORM_2, &norm);
    LOG(DEBUG) << "residual norm: " << norm;

    if (ierr)
      return false;
  }

  // dump output vector to file
  if (outputFiles)
  {
    static int evaluationNo = 0;  // counter how often this function was called

    // dumpVector(std::string filename, std::string format, Vec &vector, MPI_Comm mpiCommunicator, int componentNo=0, int nComponents=1);
    std::stringstream filename;
    filename << "out/F" << std::setw(3) << std::setfill('0') << evaluationNo++ << "r" << DihuContext::nRanksCommWorld();
    combinedVecResidual_->dumpGlobalNatural(filename.str());
  }
  assert(combinedVecResidual_->currentRepresentation() == Partition::values_representation_t::representationCombinedGlobal);
  assert(combinedVecSolution_->currentRepresentation() == Partition::values_representation_t::representationCombinedGlobal);

  // computation succeeded
  return true;
}

template<typename Term,typename MeshType, int nDisplacementComponents>
void HyperelasticitySolver<Term,MeshType,nDisplacementComponents>::
materialComputeExternalVirtualWorkDead()
{
  // compute δW_ext,dead = int_Ω B^L * phi^L * phi^M * δu^M dx + int_∂Ω T^L * phi^L * phi^M * δu^M dS

  LOG(DEBUG) << "materialComputeExternalVirtualWorkDead";

  combinedVecExternalVirtualWorkDead_->setRepresentationGlobal();
  combinedVecExternalVirtualWorkDead_->zeroEntries();               // clear entries to 0
  combinedVecExternalVirtualWorkDead_->startGhostManipulation();    // fill ghost buffers

  // get traction directly from neumann boundary conditions, they are defined element wise and the integration takes place in neumann_boundary_conditions.tpp
  std::vector<double> values;
  for (int componentNo = 0; componentNo < 3; componentNo++)
  {
    values.clear();
    neumannBoundaryConditions_->rhs()->getValuesWithoutGhosts(componentNo, values);
    LOG(DEBUG) << "component " << componentNo << ", neumannBoundaryConditions_ rhs values: " << values;


    combinedVecExternalVirtualWorkDead_->setValues(componentNo, displacementsFunctionSpace_->meshPartition()->nDofsLocalWithoutGhosts(),
                                                   displacementsFunctionSpace_->meshPartition()->dofNosLocal().data(), values.data(), INSERT_VALUES);
  }

  // integrate to account for body forces
  // -------------------------------------
  if (fabs(constantBodyForce_[0]) > 1e-12 || fabs(constantBodyForce_[1]) > 1e-12 || fabs(constantBodyForce_[2]) > 1e-12)
  {
    LOG(DEBUG) << "add contribution of body force " << constantBodyForce_;

    const int D = 3;  // dimension
    std::shared_ptr<DisplacementsFunctionSpace> functionSpace = this->data_.displacementsFunctionSpace();
    const int nDisplacementsDofsPerElement = DisplacementsFunctionSpace::nDofsPerElement();
    const int nUnknowsPerElement = nDisplacementsDofsPerElement*D;    // D directions for displacements per dof
    const int nDofsPerElement = DisplacementsFunctionSpace::nDofsPerElement();

    // define shortcuts for quadrature
    typedef Quadrature::TensorProduct<D,Quadrature::Gauss<3>> QuadratureDD;   // quadratic*quadratic = 4th order polynomial, 3 gauss points = 2*3-1 = 5th order exact

    // define type to hold evaluations of integrand
    typedef std::array<double, 3*nUnknowsPerElement> EvaluationsType;
    std::array<EvaluationsType, QuadratureDD::numberEvaluations()> evaluationsArray{};

    // setup arrays used for integration
    std::array<Vec3, QuadratureDD::numberEvaluations()> samplingPoints = QuadratureDD::samplingPoints();

    // initialize variables
    functionSpace->geometryField().setRepresentationGlobal();
    functionSpace->geometryField().startGhostManipulation();   // ensure that local ghost values of geometry field are set

    // loop over elements
    for (element_no_t elementNoLocal = 0; elementNoLocal < functionSpace->nElementsLocal(); elementNoLocal++)
    {
      // get geometry field values
      std::array<Vec3,DisplacementsFunctionSpace::nDofsPerElement()> geometry;
      functionSpace->getElementGeometry(elementNoLocal, geometry);

      // In case it is required to have a varying body force, get the body force values here, also update the integrand
      //std::array<Vec3,nDisplacementsDofsPerElement> displacementsValues;
      //this->data_.bodyForce()->getElementValues(elementNoLocal, bodyForceValues);

      // evaluate integrand at sampling points
      for (unsigned int samplingPointIndex = 0; samplingPointIndex < samplingPoints.size(); samplingPointIndex++)
      {
        // evaluate function to integrate at samplingPoints[i*2], write value to evaluations[i]
        std::array<double,D> xi = samplingPoints[samplingPointIndex];

        // compute the 3xD jacobian of the parameter space to world space mapping
        auto jacobian = DisplacementsFunctionSpace::computeJacobian(geometry, xi);
        double integrationFactor = MathUtility::computeIntegrationFactor(jacobian);

        for (unsigned int elementalDofNoM = 0; elementalDofNoM < nDofsPerElement; elementalDofNoM++)   // dof index M
        {
          for (int dimensionNo = 0; dimensionNo < 3; dimensionNo++)
          {
            double integrand = 0;

            for (unsigned int elementalDofNoL = 0; elementalDofNoL < nDofsPerElement; elementalDofNoL++)   // dof index L
            {
              integrand += DisplacementsFunctionSpace::phi(elementalDofNoL,xi) * DisplacementsFunctionSpace::phi(elementalDofNoM,xi)
                * constantBodyForce_[dimensionNo];
            }

            evaluationsArray[samplingPointIndex][elementalDofNoM*3 + dimensionNo] = integrand * integrationFactor;
          }
        }
      }  // function evaluations

      // integrate all values for the (M,d) dofs at once
      EvaluationsType integratedValues = QuadratureDD::computeIntegral(evaluationsArray);

      std::array<dof_no_t,nDisplacementsDofsPerElement> dofNosLocal = functionSpace->getElementDofNosLocal(elementNoLocal);

      // add integrated entries to result vector
      for (unsigned int elementalDofNoM = 0; elementalDofNoM < nDofsPerElement; elementalDofNoM++)   // dof index M
      {
        for (int dimensionNo = 0; dimensionNo < 3; dimensionNo++)
        {

          // get local dof no, elementalDofNoM is the dof within the element, dofNoLocal is the dof within the local subdomain
          dof_no_t dofNoLocal = dofNosLocal[elementalDofNoM];

          double integratedValue = integratedValues[elementalDofNoM*3 + dimensionNo];
          combinedVecExternalVirtualWorkDead_->setValue(dimensionNo, dofNoLocal, integratedValue, ADD_VALUES);

        }  // D
      }  // L
    }  // elementNoLocal
  }

  combinedVecExternalVirtualWorkDead_->finishGhostManipulation();     // communicate and add up values in ghost buffers
  combinedVecExternalVirtualWorkDead_->startGhostManipulation();      // communicate ghost buffers values back in place

  if (combinedVecExternalVirtualWorkDead_->containsNanOrInf())
  {
    LOG(FATAL) << "The external virtual work, δW_ext,dead contains nan or inf values. " << std::endl
      << "Check that the constantBodyForce (" << constantBodyForce_ << ") and the Neumann boundary condition values are valid.";
  }

  LOG(DEBUG) << "combinedVecExternalVirtualWorkDead (components 3-6 should be empty): " << combinedVecExternalVirtualWorkDead_->getString();
  //combinedVecExternalVirtualWorkDead_->startGhostManipulation();
}

template<typename Term,typename MeshType, int nDisplacementComponents>
void HyperelasticitySolver<Term,MeshType,nDisplacementComponents>::
materialAddAccelerationTermAndVelocityEquation(bool communicateGhosts)
{
  assert (nDisplacementComponents == 6);

  // add to solverVariableResidual_ (=combinedVecResidual_) += int_Ω rho_0 * (v^(n+1),L - v^(n),L) / dt * phi^L * phi^M * δu^M dx
  // solverVariableSolution_

  //combinedVecSolution_->startGhostManipulation();
  if (communicateGhosts)
  {
    combinedVecResidual_->startGhostManipulation();      // communicate ghost buffers values back in place
    combinedVecResidual_->zeroGhostBuffer();      // communicate ghost buffers values back in place
  }

  const int D = 3;  // dimension
  std::shared_ptr<DisplacementsFunctionSpace> functionSpace = this->data_.displacementsFunctionSpace();
  const int nDisplacementsDofsPerElement = DisplacementsFunctionSpace::nDofsPerElement();
  const int nUnknowsPerElement = nDisplacementsDofsPerElement*D;    // D directions for displacements per dof
  const int nDofsPerElement = DisplacementsFunctionSpace::nDofsPerElement();

  // define shortcuts for quadrature
  typedef Quadrature::TensorProduct<D,Quadrature::Gauss<3>> QuadratureDD;   // quadratic*quadratic = 4th order polynomial, 3 gauss points = 2*3-1 = 5th order exact

  // define type to hold evaluations of integrand
  typedef std::array<double_v_t, 3*nUnknowsPerElement> EvaluationsType;
  std::array<EvaluationsType, QuadratureDD::numberEvaluations()> evaluationsArray{};

  // setup arrays used for integration
  std::array<Vec3, QuadratureDD::numberEvaluations()> samplingPoints = QuadratureDD::samplingPoints();

  const int nElementsLocal = functionSpace->nElementsLocal();

  // loop over elements, always 4 elements at once using the vectorized functions
  for (int elementNoLocal = 0; elementNoLocal < nElementsLocal; elementNoLocal += nVcComponents)
  {

#ifdef USE_VECTORIZED_FE_MATRIX_ASSEMBLY
    // get indices of elementNos that should be handled in the current iterations,
    // this is, e.g.
    //    [10,11,12,13,-1,-1,-1,-1] (if nVcComponents==4 and nElementsLocal > 13)
    // or [10,11,12,-1,-1,-1,-1,-1] (if nVcComponents==4 and nElementsLocal == 13)

    dof_no_v_t elementNoLocalv([elementNoLocal, nElementsLocal](int i)
    {
      return (i >= nVcComponents || elementNoLocal+i >= nElementsLocal? -1: elementNoLocal+i);
    });

    // here, elementNoLocalv is the list of indices of the current iteration, e.g. [10,11,12,13,-1,-1,-1,-1]
    // elementNoLocal is the first entry of elementNoLocalv
#else
    int elementNoLocalv = elementNoLocal;
#endif

    std::array<dof_no_v_t,nDisplacementsDofsPerElement> dofNosLocal = functionSpace->getElementDofNosLocal(elementNoLocalv);

    // get geometry field values
    std::array<Vec3_v_t,DisplacementsFunctionSpace::nDofsPerElement()> geometry;
    functionSpace->getElementGeometry(elementNoLocalv, geometry);

    // get the old displacements values at the previous timestep for all dofs of the element
    std::array<Vec3_v_t,nDisplacementsDofsPerElement> oldDisplacementValues;
    this->data_.displacementsPreviousTimestep()->getElementValues(elementNoLocalv, oldDisplacementValues);

    // get the old velocity values at the previous timestep  for all dofs of the element
    std::array<Vec3_v_t,nDisplacementsDofsPerElement> oldVelocityValues;
    this->data_.velocitiesPreviousTimestep()->getElementValues(elementNoLocalv, oldVelocityValues);

    // get the old displacements values at the previous timestep for all dofs of the element
    std::array<Vec3_v_t,nDisplacementsDofsPerElement> newDisplacementValues;
    this->data_.displacements()->getElementValues(elementNoLocalv, newDisplacementValues);

    // get the old velocity values at the previous timestep  for all dofs of the element
    std::array<Vec3_v_t,nDisplacementsDofsPerElement> newVelocityValues;
    this->data_.velocities()->getElementValues(elementNoLocalv, newVelocityValues);

    // debugging check
#if 0
    std::vector<dof_no_v_t> dofNosLocalWithoutGhosts;
    functionSpace->getElementDofNosLocalWithoutGhosts(elementNoLocal, dofNosLocalWithoutGhosts);

    //if reproducing this debugging check, add calls to startGhostManipulation and finishGhostManipulation on combinedVecSolution_
    // get the new displacement values for all dofs of the element
    std::array<std::array<double,nDisplacementsDofsPerElement>,3> newDisplacementValues1;
    combinedVecSolution_->getValues(0, dofNosLocalWithoutGhosts.size(), dofNosLocalWithoutGhosts.data(), newDisplacementValues1[0].data());
    combinedVecSolution_->getValues(1, dofNosLocalWithoutGhosts.size(), dofNosLocalWithoutGhosts.data(), newDisplacementValues1[1].data());
    combinedVecSolution_->getValues(2, dofNosLocalWithoutGhosts.size(), dofNosLocalWithoutGhosts.data(), newDisplacementValues1[2].data());

    // get the new velocity values for all dofs of the element
    std::array<std::array<double,nDisplacementsDofsPerElement>,3> newVelocityValues1;
    combinedVecSolution_->getValues(3, dofNosLocalWithoutGhosts.size(), dofNosLocalWithoutGhosts.data(), newVelocityValues1[0].data());
    combinedVecSolution_->getValues(4, dofNosLocalWithoutGhosts.size(), dofNosLocalWithoutGhosts.data(), newVelocityValues1[1].data());
    combinedVecSolution_->getValues(5, dofNosLocalWithoutGhosts.size(), dofNosLocalWithoutGhosts.data(), newVelocityValues1[2].data());

    // get the new displacements values, to compare with newDisplacementValues (should be the same)
    std::array<Vec3,nDisplacementsDofsPerElement> newDisplacementValues2;
    this->data_.displacements()->getElementValues(elementNoLocal, newDisplacementValues2);

    // get the new velocity values, to compare with newVelocityValues (should be the same)
    std::array<Vec3,nDisplacementsDofsPerElement> newVelocityValues2;
    this->data_.velocities()->getElementValues(elementNoLocal, newVelocityValues2);

    LOG(DEBUG) << " element " << elementNoLocal << ":";
    for (int dofIndex = 0; dofIndex < nDisplacementsDofsPerElement; dofIndex++)
    {
      dof_no_t dofNo = dofNosLocal[dofIndex];
      LOG(DEBUG) << " dof " << dofNo << ": displacements ("
        << newDisplacementValues1[0][dofIndex] << "=" << newDisplacementValues2[dofIndex][0] << ", "
        << newDisplacementValues1[1][dofIndex] << "=" << newDisplacementValues2[dofIndex][1] << ", "
        << newDisplacementValues1[2][dofIndex] << "=" << newDisplacementValues2[dofIndex][2] << "), velocities ("
        << newVelocityValues1[0][dofIndex] << "=" << newVelocityValues2[dofIndex][0] << ", "
        << newVelocityValues1[1][dofIndex] << "=" << newVelocityValues2[dofIndex][1] << ", "
        << newVelocityValues1[2][dofIndex] << "=" << newVelocityValues2[dofIndex][2] << ")";
    }
#endif

    // evaluate integrand at sampling points
    for (unsigned int samplingPointIndex = 0; samplingPointIndex < samplingPoints.size(); samplingPointIndex++)
    {
      // evaluate function to integrate at samplingPoints[i*2], write value to evaluations[i]
      std::array<double,D> xi = samplingPoints[samplingPointIndex];

      // compute the 3xD jacobian of the parameter space to world space mapping
      auto jacobian = DisplacementsFunctionSpace::computeJacobian(geometry, xi);
      double_v_t integrationFactor = MathUtility::computeIntegrationFactor(jacobian);

      // loop over elemantal dofs, M
      for (unsigned int elementalDofNoM = 0; elementalDofNoM < nDofsPerElement; elementalDofNoM++)   // dof index M
      {
        for (int dimensionNo = 0; dimensionNo < 3; dimensionNo++)
        {
          double_v_t integrand = 0;

          // loop over elemental dofs, L
          for (unsigned int elementalDofNoL = 0; elementalDofNoL < nDofsPerElement; elementalDofNoL++)   // dof index L
          {

            const double_v_t oldVelocity = oldVelocityValues[elementalDofNoL][dimensionNo];
            const double_v_t newVelocity = newVelocityValues[elementalDofNoL][dimensionNo];

            integrand += density_ * (newVelocity - oldVelocity) / this->timeStepWidth_ * DisplacementsFunctionSpace::phi(elementalDofNoL,xi)
              * DisplacementsFunctionSpace::phi(elementalDofNoM,xi);

            evaluationsArray[samplingPointIndex][elementalDofNoM*3 + dimensionNo] = integrand * integrationFactor;
          }

          evaluationsArray[samplingPointIndex][elementalDofNoM*3 + dimensionNo] = integrand * integrationFactor;
        }
      }
    }  // function evaluations

    // integrate all values for the (M,d) dofs at once
    EvaluationsType integratedValues = QuadratureDD::computeIntegral(evaluationsArray);

    for (unsigned int elementalDofNoM = 0; elementalDofNoM < nDofsPerElement; elementalDofNoM++)   // dof index M
    {
      for (int dimensionNo = 0; dimensionNo < 3; dimensionNo++)
      {

        // get local dof no, elementalDofNoM is the dof within the element, dofNoLocal is the dof within the local subdomain
        dof_no_v_t dofNoLocal = dofNosLocal[elementalDofNoM];

        // add integrated entries of velocityTerm to result vector
        double_v_t integratedValue = integratedValues[elementalDofNoM*3 + dimensionNo];
        combinedVecResidual_->setValue(dimensionNo, dofNoLocal, integratedValue, ADD_VALUES);

#if 0
        if (VLOG_IS_ON(1))
        {
          global_no_t dofNoGlobal = functionSpace->meshPartition()->getDofNoGlobalPetsc(dofNoLocal);
          global_no_t elementNoGlobal = functionSpace->meshPartition()->getElementNoGlobalNatural(elementNoLocal);
          VLOG(1) << "el global " << elementNoGlobal << " add to dof global " << dofNoGlobal << " value " << integratedValue;
        }
#endif
        // compute velocity equation
        const double_v_t oldDisplacement = oldDisplacementValues[elementalDofNoM][dimensionNo];
        const double_v_t newDisplacement = newDisplacementValues[elementalDofNoM][dimensionNo];

        //double_v_t oldVelocity = oldVelocityValues[elementalDofNoM][dimensionNo];
        const double_v_t newVelocity = newVelocityValues[elementalDofNoM][dimensionNo];

        // add the velocity/displacement equation 1/dt (u^(n+1) - u^(n)) - v^(n+1) = 0 in the velocity slot
        double_v_t residuum = 1.0 / this->timeStepWidth_ * (newDisplacement - oldDisplacement) - newVelocity;

        // only set value if current dof is local, we must not set dof values here, even when "INSERT_VALUES" is used, because the dof values would still be summed up by finishGhostManipulation() (even though it is INSERT_VALUES and not ADD_VALUES)
#ifdef USE_VECTORIZED_FE_MATRIX_ASSEMBLY
        for (int vcComponentNo = 0; vcComponentNo < Vc::double_v::size(); vcComponentNo++)
        {
          if (dofNoLocal[vcComponentNo] != -1 && dofNoLocal[vcComponentNo] < functionSpace->nDofsLocalWithoutGhosts())
          {
            combinedVecResidual_->setValue(3 + dimensionNo, dofNoLocal[vcComponentNo], residuum[vcComponentNo], INSERT_VALUES);
          }
        }
#else
        if (dofNoLocal < functionSpace->nDofsLocalWithoutGhosts())
        {
          combinedVecResidual_->setValue(3 + dimensionNo, dofNoLocal, residuum, INSERT_VALUES);
        }
#endif

      }  // D
    }  // L
  }  // elementNoLocal

  if (communicateGhosts)
  {
    combinedVecResidual_->finishGhostManipulation();     // communicate and add up values in ghost buffers
  }

  //combinedVecSolution_->zeroGhostBuffer();
  //combinedVecSolution_->finishGhostManipulation();
  //LOG(DEBUG) << "combinedVecExternalVirtualWorkDead: " << combinedVecExternalVirtualWorkDead_->getString();

  //combinedVecExternalVirtualWorkDead_->startGhostManipulation();
}

template<typename Term,typename MeshType, int nDisplacementComponents>
bool HyperelasticitySolver<Term,MeshType,nDisplacementComponents>::
materialComputeJacobian()
{
  // analytic jacobian combinedMatrixJacobian_
  //  output is combinedMatrixJacobian_, a PartitionedMatHyperelasticity or solverMatrixJacobian_, the normal Mat, contains no Dirichlet BC dofs
  //  input is solverVariableSolution_, a normal Vec, the same values have already been assigned to this->data_.displacements() and this->data_.pressure()

  const bool outputValues = false;
  if (outputValues)
    LOG(DEBUG) << "input: " << getString(solverVariableSolution_);

  // assert that data representation is global
  assert(combinedVecSolution_->currentRepresentation() == Partition::values_representation_t::representationCombinedGlobal);

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

  // define types to hold evaluations of integrand
  typedef std::array<double_v_t, nUnknowsPerElement*nUnknowsPerElement> EvaluationsDisplacementsType;
  std::array<EvaluationsDisplacementsType, QuadratureDD::numberEvaluations()> evaluationsArrayDisplacements{};

  typedef std::array<double_v_t, nPressureDofsPerElement*nUnknowsPerElement> EvaluationsPressureType;
  std::array<EvaluationsPressureType, QuadratureDD::numberEvaluations()> evaluationsArrayPressure{};

  typedef std::array<double_v_t, nDisplacementsDofsPerElement*nDisplacementsDofsPerElement> EvaluationsUVType;
  std::array<EvaluationsUVType, QuadratureDD::numberEvaluations()> evaluationsArrayUV{};

  // setup arrays used for integration
  std::array<Vec3, QuadratureDD::numberEvaluations()> samplingPoints = QuadratureDD::samplingPoints();

  // loop over elements, always 4 elements at once using the vectorized functions
  for (int elementNoLocal = 0; elementNoLocal < nElementsLocal; elementNoLocal += nVcComponents)
  {

#ifdef USE_VECTORIZED_FE_MATRIX_ASSEMBLY
    // get indices of elementNos that should be handled in the current iterations,
    // this is, e.g.
    //    [10,11,12,13,-1,-1,-1,-1] (if nVcComponents==4 and nElementsLocal > 13)
    // or [10,11,12,-1,-1,-1,-1,-1] (if nVcComponents==4 and nElementsLocal == 13)

    dof_no_v_t elementNoLocalv([elementNoLocal, nElementsLocal](int i)
    {
      return (i >= nVcComponents || elementNoLocal+i >= nElementsLocal? -1: elementNoLocal+i);
    });

    // here, elementNoLocalv is the list of indices of the current iteration, e.g. [10,11,12,13,-1,-1,-1,-1]
    // elementNoLocal is the first entry of elementNoLocalv
#else
    int elementNoLocalv = elementNoLocal;
#endif

    // get indices of element-local dofs
    std::array<dof_no_v_t,nDisplacementsDofsPerElement> dofNosLocal = displacementsFunctionSpace->getElementDofNosLocal(elementNoLocalv);
    std::array<dof_no_v_t,nPressureDofsPerElement> dofNosLocalPressure = pressureFunctionSpace->getElementDofNosLocal(elementNoLocalv);

    // set values to zero to be able to add values later

    // initialize top-left submatrix to zero, for dynamic case initialize submatrices (1,0), (0,1), (1,1)
    // loop over indices of unknows (aDof,aComponent),(bDof,bComponent) or (L,a),(M,b)
    for (int aDof = 0; aDof < nDisplacementsDofsPerElement; aDof++)        // L
    {
      for (int aComponent = 0; aComponent < D; aComponent++)               // a
      {
        for (int bDof = 0; bDof < nDisplacementsDofsPerElement; bDof++)    // M
        {
          for (int bComponent = 0; bComponent < D; bComponent++)           // b
          {
            // get local dof no, aDof is the dof within the element, dofNoLocal is the dof within the local subdomain
            dof_no_v_t dofANoLocal = dofNosLocal[aDof];
            dof_no_v_t dofBNoLocal = dofNosLocal[bDof];

            combinedMatrixJacobian_->setValue(aComponent, dofANoLocal, bComponent, dofBNoLocal, 0.0, INSERT_VALUES);

            // for dynamic case also initialize center-left, top-center and center-center sub matrices
            if (nDisplacementComponents == 6)
            {
              // set entry of top-center (0,1) sub matrix, l_δu,Δv
              // this entry will be computed by an integral
              combinedMatrixJacobian_->setValue(aComponent, dofANoLocal, 3+bComponent, dofBNoLocal, 0.0, INSERT_VALUES);

              // set entry of center-left (1,0) sub matrix, l_δv,Δu
              // this entry can directly be computed
              const int delta_ab = (aComponent == bComponent? 1 : 0);
              const int delta_LM = (aDof == bDof? 1 : 0);

              const double entryVU = 1./this->timeStepWidth_ * delta_ab * delta_LM;
              combinedMatrixJacobian_->setValue(3+aComponent, dofANoLocal, bComponent, dofBNoLocal, entryVU, INSERT_VALUES);

              // set entry of center-center (1,1) sub matrix, l_δv,Δv
              // this entry can directly be computed
              const double entryVV = -delta_ab * delta_LM;
              combinedMatrixJacobian_->setValue(3+aComponent, dofANoLocal, 3+bComponent, dofBNoLocal, entryVV, INSERT_VALUES);
            }
          }  // b
        }  // M
      }  // a
    }  // L

    // initialize top right and bottom left sub matrices to zero, those entries are only present for the incompressible formulation
    if (Term::isIncompressible)
    {
      // loop over indices of unknows aDof,(bDof,bComponent) or L,(M,b)
      for (int lDof = 0; lDof < nPressureDofsPerElement; lDof++)           // L
      {
        for (int aDof = 0; aDof < nDisplacementsDofsPerElement; aDof++)    // M
        {
          for (int aComponent = 0; aComponent < D; aComponent++)           // a
          {
            // get local dof no, aDof is the dof within the element, dofNoLocal is the dof within the local subdomain
            dof_no_v_t dofLNoLocal = dofNosLocalPressure[lDof];     // dof with respect to pressure function space
            dof_no_v_t dofMNoLocal = dofNosLocal[aDof];             // dof with respect to displacements function space

            const int pressureDofNo = nDisplacementComponents;  // 3 or 6, depending if static or dynamic problem

            // set entry in lower left submatrix
            // parameters: componentNoRow, dofNoLocalRow, componentNoColumn, dofNoLocalColumn, value
            combinedMatrixJacobian_->setValue(pressureDofNo, dofLNoLocal, aComponent, dofMNoLocal, 0.0, INSERT_VALUES);

            // set entry in upper right submatrix
            combinedMatrixJacobian_->setValue(aComponent, dofMNoLocal, pressureDofNo, dofLNoLocal, 0.0, INSERT_VALUES);

          }  // aComponent
        }  // aDof
      }  // lDof

      // loop over diagonal matrix entries in p-part (bottom right submatrix), set diagonal entries to 0
      // This allocates nonzero entries and sets them to zero. It is needed for the solver.

      // The direct solvers are not able to solve the saddle point problem in serial execution, add an artifical regularization term.
      double epsilon = 1e-12;

      // or parallel execution everything is fine and we do not need regularization
      if (this->data_.functionSpace()->meshPartition()->nRanks() > 1)
        epsilon = 0;

      for (int lDof = 0; lDof < nPressureDofsPerElement; lDof++)           // L
      {
        const int pressureDofNo = nDisplacementComponents;  // 3 or 6, depending if static or dynamic problem

        dof_no_v_t dofLNoLocal = dofNosLocalPressure[lDof];     // dof with respect to pressure function space
        combinedMatrixJacobian_->setValue(pressureDofNo, dofLNoLocal, pressureDofNo, dofLNoLocal, epsilon, INSERT_VALUES);
      }
    }  // elementNoLocal

    // set diagonal to zero, this would be correct but for some reason the solvers do not like systems with zero diagonal, therefore epsilon was set on the diagonal
    //PetscErrorCode ierr;
    //ierr = MatDiagonalSet(combinedMatrixJacobian_->valuesGlobal(), zeros_, INSERT_VALUES); CHKERRQ(ierr);

  }  // if Term::isIncompressible

  // allow switching between stiffnessMatrix->setValue(... INSERT_VALUES) and ADD_VALUES
  combinedMatrixJacobian_->assembly(MAT_FLUSH_ASSEMBLY);

  // loop over elements, always 4 elements at once using the vectorized functions
  for (int elementNoLocal = 0; elementNoLocal < nElementsLocal; elementNoLocal += nVcComponents)
  {

#ifdef USE_VECTORIZED_FE_MATRIX_ASSEMBLY
    // get indices of elementNos that should be handled in the current iterations,
    // this is, e.g.
    //    [10,11,12,13,-1,-1,-1,-1] (if nVcComponents==4 and nElementsLocal > 13)
    // or [10,11,12,-1,-1,-1,-1,-1] (if nVcComponents==4 and nElementsLocal == 13)

    dof_no_v_t elementNoLocalv([elementNoLocal, nElementsLocal](dof_no_t i)
    {
      return (i >= nVcComponents || elementNoLocal+i >= nElementsLocal? -1: elementNoLocal+i);
    });

    // here, elementNoLocalv is the list of indices of the current iteration, e.g. [10,11,12,13,-1,-1,-1,-1]
    // elementNoLocal is the first entry of elementNoLocalv
#else
    int elementNoLocalv = elementNoLocal;
#endif

    // get geometry field of reference configuration
    std::array<Vec3_v_t,nDisplacementsDofsPerElement> geometryReferenceValues;
    this->data_.geometryReference()->getElementValues(elementNoLocalv, geometryReferenceValues);

    // get displacements field values for element
    std::array<Vec3_v_t,nDisplacementsDofsPerElement> displacementsValues;
    this->data_.displacements()->getElementValues(elementNoLocalv, displacementsValues);

    //LOG(DEBUG) << "elementNoLocal " << elementNoLocal << ", displacementsValues: " << displacementsValues;

    std::array<double_v_t,nPressureDofsPerElement> pressureValuesCurrentElement;
    if (Term::isIncompressible)
      this->data_.pressure()->getElementValues(elementNoLocalv, pressureValuesCurrentElement);

    std::array<Vec3_v_t,nDisplacementsDofsPerElement> elementalDirectionValues;
    this->data_.fiberDirection()->getElementValues(elementNoLocalv, elementalDirectionValues);

    // loop over integration points (e.g. gauss points) for displacements field
    for (unsigned int samplingPointIndex = 0; samplingPointIndex < samplingPoints.size(); samplingPointIndex++)
    {
      // get parameter values of current sampling point
      Vec3 xi = samplingPoints[samplingPointIndex];

      // compute the 3x3 jacobian of the parameter space to world space mapping
      Tensor2_v_t<D> jacobianMaterial = DisplacementsFunctionSpace::computeJacobian(geometryReferenceValues, xi);
      double_v_t jacobianDeterminant;
      Tensor2_v_t<D> inverseJacobianMaterial = MathUtility::computeInverse(jacobianMaterial, jacobianDeterminant);

      // jacobianMaterial[columnIdx][rowIdx] = dX_rowIdx/dxi_columnIdx
      // inverseJacobianMaterial[columnIdx][rowIdx] = dxi_rowIdx/dX_columnIdx because of inverse function theorem

      // get the factor in the integral that arises from the change in integration domain from world to parameter space
      double_v_t integrationFactor = MathUtility::abs(jacobianDeterminant);   //MathUtility::computeIntegrationFactor(jacobianMaterial);

      Tensor2_v_t<D> deformationGradient = this->computeDeformationGradient(displacementsValues, inverseJacobianMaterial, xi);    // F
      double_v_t deformationGradientDeterminant;    // J
      Tensor2_v_t<D> inverseDeformationGradient = MathUtility::computeInverse(deformationGradient, deformationGradientDeterminant);  // F^-1

      Tensor2_v_t<D> rightCauchyGreen = this->computeRightCauchyGreenTensor(deformationGradient);  // C = F^T*F

      double_v_t rightCauchyGreenDeterminant;   // J^2
      Tensor2_v_t<D> inverseRightCauchyGreen = MathUtility::computeSymmetricInverse(rightCauchyGreen, rightCauchyGreenDeterminant);  // C^-1

      // fiber direction
      Vec3_v_t fiberDirection = displacementsFunctionSpace->template interpolateValueInElement<3>(elementalDirectionValues, xi);

#ifndef NDEBUG
      if (Term::usesFiberDirection)
      {
        if (fabs(MathUtility::norm<3>(fiberDirection) - 1) > 1e-3)
          LOG(FATAL) << "fiberDirecton " << fiberDirection << " is not normalized, elementalDirectionValues:" << elementalDirectionValues;
      }
#endif

      // invariants
      std::array<double_v_t,5> invariants = this->computeInvariants(rightCauchyGreen, rightCauchyGreenDeterminant, fiberDirection);  // I_1, I_2, I_3
      std::array<double_v_t,5> reducedInvariants = this->computeReducedInvariants(invariants, deformationGradientDeterminant); // Ibar_1, ..., Ibar_5

      // pressure is the separately interpolated pressure for mixed formulation
      double_v_t pressure = 0;
      if (Term::isIncompressible)
        pressure = pressureFunctionSpace->interpolateValueInElement(pressureValuesCurrentElement, xi);

      // Pk2 stress tensor S = S_vol + S_iso (p.234)
      //! compute 2nd Piola-Kirchhoff stress tensor S = 2*dPsi/dC and the fictitious PK2 Stress Sbar
      Tensor2_v_t<D> fictitiousPK2Stress;   // Sbar
      Tensor2_v_t<D> pk2StressIsochoric;    // S_iso
      Tensor2_v_t<D> pK2Stress = this->computePK2Stress(pressure, rightCauchyGreen, inverseRightCauchyGreen, invariants, reducedInvariants,
                                                        deformationGradientDeterminant, fiberDirection,
                                                        fictitiousPK2Stress, pk2StressIsochoric);

      std::array<Vec3,nDisplacementsDofsPerElement> gradPhi = displacementsFunctionSpace->getGradPhi(xi);
      // (column-major storage) gradPhi[L][a] = dphi_L / dxi_a
      // gradPhi[column][row] = gradPhi[dofIndex][i] = dphi_dofIndex/dxi_i, columnIdx = dofIndex, rowIdx = which direction


      Tensor4_v_t<D> elasticityTensor;
      Tensor4_v_t<D> fictitiousElasticityTensor;
      Tensor4_v_t<3> elasticityTensorIso;
      computeElasticityTensor(rightCauchyGreen, inverseRightCauchyGreen, deformationGradientDeterminant, pressure, invariants, reducedInvariants, fictitiousPK2Stress, pk2StressIsochoric, fiberDirection,
                              fictitiousElasticityTensor, elasticityTensorIso, elasticityTensor);

      // test if implementation of S is correct
      this->materialTesting(pressure, rightCauchyGreen, inverseRightCauchyGreen, reducedInvariants, deformationGradientDeterminant, fiberDirection, fictitiousPK2Stress, pk2StressIsochoric);

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
      VLOG(2) << "  pK2Stress: S=" << pK2Stress;
      VLOG(2) << "  gradPhi: " << gradPhi;

      VLOG(1) << "  sampling point " << samplingPointIndex << "/" << samplingPoints.size() << ", xi: " << xi << ", J: " << deformationGradientDeterminant << ", p: " << pressure << ", S11: " << pK2Stress[0][0];

      if (Vc::any_of(deformationGradientDeterminant < 1e-12))   // if any entry of the deformation gradient is negative
      {
        LOG(WARNING) << "Deformation gradient " << deformationGradient << " has zero or negative determinant " << deformationGradientDeterminant
          << std::endl << "Geometry values in element " << elementNoLocal << ": " << geometryReferenceValues << std::endl
          << "Displacements at xi " << xi << ": " << displacementsValues;

        lastSolveSucceeded_ = false;
      }

      // add contributions of submatrix uu (upper left)

      // loop over pairs basis functions and evaluate integrand at xi
      for (int aDof = 0; aDof < nDisplacementsDofsPerElement; aDof++)    // index over dofs, each dof has D components, L in derivation
      {
        for (int aComponent = 0; aComponent < D; aComponent++)     // lower-case a in derivation, index over displacements components
        {

          for (int bDof = 0; bDof < nDisplacementsDofsPerElement; bDof++)  // index over dofs, each dof has D components, M in derivation
          {
            for (int bComponent = 0; bComponent < D; bComponent++)     // lower-case b in derivation, index over displacements components
            {
              double_v_t integrand = 0.0;

              for (int bInternal = 0; bInternal < D; bInternal++)     // capital B in derivation
              {
                for (int dInternal = 0; dInternal < D; dInternal++)     // capital D in derivation
                {
                  // compute integrand phi_La,B * tilde{k}_abBD * phi_Mb,D

                  // ----------------------------
                  // compute derivatives of phi
                  double_v_t dphiL_dXB = 0.0;
                  double_v_t dphiM_dXD = 0.0;

                  // helper index k for multiplication with inverse Jacobian
                  for (int k = 0; k < D; k++)
                  {
                    // (column-major storage) gradPhi[L][k] = dphi_L / dxi_k
                    // gradPhi[column][row] = gradPhi[dofIndex][k] = dphi_dofIndex/dxi_k, columnIdx = dofIndex, rowIdx = which direction

                    // compute dphiL/dXB from dphiL/dxik and dxik/dXB
                    const double dphiL_dxik = gradPhi[aDof][k];    // dphi_L/dxik
                    const double_v_t dxik_dXB = inverseJacobianMaterial[bInternal][k];  // inverseJacobianMaterial[B][k] = J^{-1}_kB = dxi_k/dX_B

                    dphiL_dXB += dphiL_dxik * dxik_dXB;

                    // compute dphiM/dXD from dphiM/dxik and dxik/dXD
                    const double dphiM_dxik = gradPhi[bDof][k];    // dphi_M/dxik
                    const double_v_t dxik_dXD = inverseJacobianMaterial[dInternal][k];  // inverseJacobianMaterial[D][k] = J^{-1}_kD = dxi_k/dX_D

                    dphiM_dXD += dphiM_dxik * dxik_dXD;
                  }   // k

                  const double_v_t sBD = pK2Stress[dInternal][bInternal];
                  const int delta_ab = (aComponent == bComponent? 1 : 0);

                  double_v_t k_abBD = delta_ab * sBD;

                  for (int cInternal = 0; cInternal < D; cInternal++)     // capital C in derivation
                  {
                    for (int aInternal = 0; aInternal < D; aInternal++)     // capital A in derivation
                    {
                      const double_v_t faA = deformationGradient[aInternal][aComponent];
                      const double_v_t fbC = deformationGradient[cInternal][bComponent];

                      const double_v_t cABCD = elasticityTensor[dInternal][cInternal][bInternal][aInternal];  // get c_{ABCD}

                      k_abBD += faA * fbC * cABCD;
                    }   // A
                  }   // C

                  integrand += dphiL_dXB * k_abBD * dphiM_dXD;

                }  // D
              }  // B

              VLOG(2) << "   (L,a)=(" << aDof << "," << aComponent << "), integrand: " << integrand;

              // compute index of degree of freedom and component (result vector index)
              const int j = aDof*D + aComponent;
              const int i = bDof*D + bComponent;
              const int index = j*nUnknowsPerElement + i;

              // store integrand in evaluations array
              evaluationsArrayDisplacements[samplingPointIndex][index] = integrand * integrationFactor;

            }  // b, bComponent
          }   // M, bDof
        }  // a, aComponent
      }  // L, aDof

      // add contributions of submatrix up and pu (lower left and upper right), only for incompressible formulation
      if (Term::isIncompressible)
      {
        // loop over indices of unknows aDof,(bDof,bComponent) or L,(M,b)
        for (int lDof = 0; lDof < nPressureDofsPerElement; lDof++)           // L
        {
          for (int aDof = 0; aDof < nDisplacementsDofsPerElement; aDof++)    // M
          {
            for (int aComponent = 0; aComponent < D; aComponent++)           // a
            {
              double_v_t fInv_Ba_dphiM_dXB = 0.0;

              for (int bInternal = 0; bInternal < D; bInternal++)     // capital B in derivation
              {
                // compute derivatives of phi
                double_v_t dphiM_dXB = 0.0;

                // helper index k for multiplication with inverse Jacobian
                for (int k = 0; k < D; k++)
                {
                  // (column-major storage) gradPhi[L][k] = dphi_L / dxi_k
                  // gradPhi[column][row] = gradPhi[dofIndex][k] = dphi_dofIndex/dxi_k, columnIdx = dofIndex, rowIdx = which direction

                  // compute dphiM/dXB from dphiM/dxik and dxik/dXB
                  const double dphiM_dxik = gradPhi[aDof][k];    // dphi_M/dxik
                  const double_v_t dxik_dXB = inverseJacobianMaterial[bInternal][k];  // inverseJacobianMaterial[B][k] = J^{-1}_kB = dxi_k/dX_B

                  dphiM_dXB += dphiM_dxik * dxik_dXB;
                }   // k

                const double_v_t fInv_Ba = inverseDeformationGradient[aComponent][bInternal];

                fInv_Ba_dphiM_dXB += fInv_Ba * dphiM_dXB;
              }

              // compute integrand J * psi_L * (F^-1)_Ba * phi_Ma,B

              const double_v_t psiL = pressureFunctionSpace->phi(lDof,xi);
              const double_v_t integrand = deformationGradientDeterminant * psiL * fInv_Ba_dphiM_dXB;

              // compute index of degree of freedom and component (result vector index)
              const int j = lDof;
              const int i = aDof*D + aComponent;
              const int index = j*nUnknowsPerElement + i;

              // store integrand in evaluations array
              evaluationsArrayPressure[samplingPointIndex][index] = integrand * integrationFactor;

            }  // a
          }  // M
        }  // L
      }  // if incompressible

      // add contributions of submatrix uv (top-center, only for dynamic problem)
      if (nDisplacementComponents == 6)
      {
        for (int lDof = 0; lDof < nDisplacementsDofsPerElement; lDof++)    // index over dofs, each dof has D components, L in derivation
        {
          for (int mDof = 0; mDof < nDisplacementsDofsPerElement; mDof++)  // index over dofs, each dof has D components, M in derivation
          {
            // integrate ∫_Ω ρ0 ϕ^L ϕ^M dV, the actual needed value is 1/dt δ_ab ∫_Ω ρ0 ϕ^L ϕ^M dV, but this will be computed later
            const double integrand = density_ * displacementsFunctionSpace->phi(lDof, xi) * displacementsFunctionSpace->phi(mDof, xi);

            // compute index of degree of freedom and component (result vector index)
            const int index = lDof*nDisplacementsDofsPerElement + mDof;

            // store integrand in evaluations array
            evaluationsArrayUV[samplingPointIndex][index] = integrand * integrationFactor;
          }   // M, mDof
        }   // L, lDof

      }  // if dynamic problem
    }   // sampling points

    // integrate all values for result vector entries at once
    EvaluationsDisplacementsType integratedValuesDisplacements = QuadratureDD::computeIntegral(evaluationsArrayDisplacements);

    EvaluationsPressureType integratedValuesPressure;
    if (Term::isIncompressible)
    {
      integratedValuesPressure = QuadratureDD::computeIntegral(evaluationsArrayPressure);
    }

    EvaluationsUVType integratedValuesUV;
    if (nDisplacementComponents == 6)
    {
      integratedValuesUV = QuadratureDD::computeIntegral(evaluationsArrayUV);
    }

    // get indices of element-local dofs
    std::array<dof_no_v_t,nDisplacementsDofsPerElement> dofNosLocal = displacementsFunctionSpace->getElementDofNosLocal(elementNoLocalv);
    std::array<dof_no_v_t,nPressureDofsPerElement> dofNosLocalPressure = pressureFunctionSpace->getElementDofNosLocal(elementNoLocalv);

    // add entries in result stiffness matrix for displacements (upper left part)

    // loop over indices of unknows (aDof,aComponent),(bDof,bComponent) or (L,a),(M,b)
    for (int aDof = 0; aDof < nDisplacementsDofsPerElement; aDof++)        // L
    {
      for (int aComponent = 0; aComponent < D; aComponent++)               // a
      {
        for (int bDof = 0; bDof < nDisplacementsDofsPerElement; bDof++)    // M
        {
          for (int bComponent = 0; bComponent < D; bComponent++)           // b
          {
            // compute index of degree of freedom and component for array of integrated values
            const int j = aDof*D + aComponent;
            const int i = bDof*D + bComponent;
            const int index = j*nUnknowsPerElement + i;

            // integrate value and set entry
            double_v_t integratedValue = integratedValuesDisplacements[index];

            // get local dof no, aDof is the dof within the element, dofNoLocal is the dof within the local subdomain
            dof_no_v_t dofANoLocal = dofNosLocal[aDof];
            dof_no_v_t dofBNoLocal = dofNosLocal[bDof];

            VLOG(1) << "  result entry (L,a)=(" <<aDof<< "," <<aComponent<< "), (M,b)=(" <<bDof<< "," <<bComponent<< ") "
              << ", dof (" << dofANoLocal << "," << dofBNoLocal << ")"
              << ", integrated value: " <<integratedValue;
            //VLOG(1) << "  jacobian[" << displacementsFunctionSpace->meshPartition()->getDofNoGlobalPetsc(dofANoLocal) << "," << aComponent << "; "
            //  << displacementsFunctionSpace->meshPartition()->getDofNoGlobalPetsc(dofBNoLocal) << "," << bComponent << "] = " << integratedValue;

            // parameters: componentNoRow, dofNoLocalRow, componentNoColumn, dofNoLocalColumn, value
            combinedMatrixJacobian_->setValue(aComponent, dofANoLocal, bComponent, dofBNoLocal, integratedValue, ADD_VALUES);

          }  // bComponent
        }  // bDof
      }  // aComponent
    }  // aDof

    // add entries in result stiffness matrix for pressure (lower left and upper right parts, up and pu, symmetric), only for incompressible formulation
    if (Term::isIncompressible)
    {
      // loop over indices of unknows aDof,(bDof,bComponent) or L,(M,b)
      for (int lDof = 0; lDof < nPressureDofsPerElement; lDof++)           // L
      {
        for (int aDof = 0; aDof < nDisplacementsDofsPerElement; aDof++)    // M
        {
          for (int aComponent = 0; aComponent < D; aComponent++)           // a
          {
            // compute index of degree of freedom and component for array of integrated values
            const int j = lDof;
            const int i = aDof*D + aComponent;
            const int index = j*nUnknowsPerElement + i;

            // get result of quadrature
            double_v_t integratedValue = integratedValuesPressure[index];

            // get local dof no, aDof is the dof within the element, dofNoLocal is the dof within the local subdomain
            dof_no_v_t dofLNoLocal = dofNosLocalPressure[lDof];     // dof with respect to pressure function space
            dof_no_v_t dofMNoLocal = dofNosLocal[aDof];  // dof with respect to displacements function space

            // set entry in lower left submatrix

            const int pressureDofNo = nDisplacementComponents;  // 3 or 6, depending if static or dynamic problem

            // parameters: componentNoRow, dofNoLocalRow, componentNoColumn, dofNoLocalColumn, value
            combinedMatrixJacobian_->setValue(pressureDofNo, dofLNoLocal, aComponent, dofMNoLocal, integratedValue, ADD_VALUES);

            // set entry in upper right submatrix
            combinedMatrixJacobian_->setValue(aComponent, dofMNoLocal, pressureDofNo, dofLNoLocal, integratedValue, ADD_VALUES);

          }  // aComponent
        }  // aDof
      }  // lDof
    }

    // add entries in resulting stiffness matrix for submatrix uv (top-center, only for dynamic problem)
    if (nDisplacementComponents == 6)
    {
      for (int lDof = 0; lDof < nDisplacementsDofsPerElement; lDof++)    // index over dofs, each dof has D components, L in derivation
      {
        for (int mDof = 0; mDof < nDisplacementsDofsPerElement; mDof++)  // index over dofs, each dof has D components, M in derivation
        {
          // get local dof no, lDof is the dof within the element, dofNoLocal is the dof within the local subdomain
          dof_no_v_t dofLNoLocal = dofNosLocal[lDof];
          dof_no_v_t dofMNoLocal = dofNosLocal[mDof];

          // compute index
          const int index = lDof*nDisplacementsDofsPerElement + mDof;

          // get result of quadrature
          const double_v_t integratedValue = integratedValuesUV[index];

          // integratedValue is only ∫_Ω ρ0 ϕ^L ϕ^M dV,
          // but we need 1/dt δ_ab ∫_Ω ρ0 ϕ^L ϕ^M dV

          for (int aComponent = 0; aComponent < D; aComponent++)           // a
          {
            for (int bComponent = 0; bComponent < D; bComponent++)           // b
            {
              if (aComponent != bComponent)
                continue;

              double_v_t resultingValue = 1./this->timeStepWidth_ * integratedValue;

              // set entrie
              // parameters: componentNoRow, dofNoLocalRow, componentNoColumn, dofNoLocalColumn, value
              combinedMatrixJacobian_->setValue(aComponent, dofMNoLocal, 3+bComponent, dofLNoLocal, resultingValue, ADD_VALUES);

            }  // bComponent
          }  // aComponent
        }  // M
      }  // L
    }
  }  // local elements

  combinedMatrixJacobian_->assembly(MAT_FINAL_ASSEMBLY);

  if (!lastSolveSucceeded_)
  {
    // return false means computation was not successful
    return false;
  }

  // computation was successful (no negative jacobian)
  return true;
}

template<typename Term,typename MeshType, int nDisplacementComponents>
unsigned int HyperelasticitySolver<Term,MeshType,nDisplacementComponents>::
materialDetermineNumberNonzerosInJacobian()
{
  unsigned int nNonZeros = 0;

  // get pointer to function space
  std::shared_ptr<DisplacementsFunctionSpace> displacementsFunctionSpace = this->data_.displacementsFunctionSpace();
  std::shared_ptr<PressureFunctionSpace> pressureFunctionSpace = this->data_.pressureFunctionSpace();

  const int D = 3;  // dimension
  const int nDisplacementsDofsPerElement = DisplacementsFunctionSpace::nDofsPerElement();
  const int nPressureDofsPerElement = PressureFunctionSpace::nDofsPerElement();
  const int nElementsLocal = displacementsFunctionSpace->nElementsLocal();

  nNonZeros = nElementsLocal * MathUtility::sqr(nDisplacementsDofsPerElement * D);

  if (nDisplacementComponents == 6)
  {
    nNonZeros += nElementsLocal * MathUtility::sqr(nDisplacementsDofsPerElement * D) * 3;
  }

  if (Term::isIncompressible)
  {
    nNonZeros += nElementsLocal * nPressureDofsPerElement * nDisplacementsDofsPerElement * D * 2;
    nNonZeros += nElementsLocal * nPressureDofsPerElement;
  }

  return nNonZeros;
}

} // namespace
