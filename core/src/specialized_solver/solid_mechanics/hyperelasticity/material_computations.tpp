#include "specialized_solver/solid_mechanics/hyperelasticity/hyperelasticity_solver.h"

#include <Python.h>  // has to be the first included header
#include <array>
#include <Vc/Vc>

#include "equation/mooney_rivlin_incompressible.h"

namespace SpatialDiscretization
{


template<typename Term,typename MeshType, int nDisplacementComponents>
void HyperelasticitySolver<Term,MeshType,nDisplacementComponents>::
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
      VLOG(1) << "elementNoGlobal " << elementNoGlobal << ", displacementsValues: " << displacementsValues;
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

      // F
      Tensor2_v_t<D> deformationGradient = this->computeDeformationGradient(displacementsValues, inverseJacobianMaterial, xi);
      double_v_t deformationGradientDeterminant = MathUtility::computeDeterminant(deformationGradient);  // J

      Tensor2_v_t<D> rightCauchyGreen = this->computeRightCauchyGreenTensor(deformationGradient);  // C = F^T*F

      double_v_t rightCauchyGreenDeterminant;   // J^2
      Tensor2_v_t<D> inverseRightCauchyGreen = MathUtility::computeSymmetricInverse(rightCauchyGreen, rightCauchyGreenDeterminant);  // C^-1

      // fiber direction
      Vec3_v_t fiberDirection = displacementsFunctionSpace->template interpolateValueInElement<3>(elementalDirectionValues, xi);

      // invariants
      std::array<double_v_t,5> invariants = this->computeInvariants(rightCauchyGreen, rightCauchyGreenDeterminant, fiberDirection);  // I_1, I_2, I_3, I_4, I_5
      std::array<double_v_t,5> reducedInvariants = this->computeReducedInvariants(invariants, deformationGradientDeterminant); // Ibar_1, Ibar_2, Ibar_4, Ibar_5

      // pressure is the separately interpolated pressure for mixed formulation
      double_v_t pressure = pressureFunctionSpace->interpolateValueInElement(pressureValuesCurrentElement, xi);

      // Pk2 stress tensor S = S_vol + S_iso (p.234)
      //! compute 2nd Piola-Kirchhoff stress tensor S = 2*dPsi/dC and the fictitious PK2 Stress Sbar
      Tensor2_v_t<D> fictitiousPK2Stress;   // Sbar
      Tensor2_v_t<D> pk2StressIsochoric;    // S_iso
      Tensor2_v_t<D> pK2Stress = this->computePK2Stress(pressure, rightCauchyGreen, inverseRightCauchyGreen, reducedInvariants, deformationGradientDeterminant, fiberDirection,
                                                       fictitiousPK2Stress, pk2StressIsochoric
                                                    );

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
      VLOG(2) << "  PK2Stress: S=" << pK2Stress;
      VLOG(2) << "  gradPhi: " << gradPhi;

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
          evaluationsArrayDisplacements[samplingPointIndex][i] = integrand;

        }  // a
      }  // L

      // loop over basis functions and evaluate integrand at xi for pressure part ((J-1)*psi)
      for (int dofIndex = 0; dofIndex < nPressureDofsPerElement; dofIndex++)           // index over dofs in element, L in derivation
      {
        const double phiL = pressureFunctionSpace->phi(dofIndex, xi);
        const double_v_t integrand = (deformationGradientDeterminant - 1.0) * phiL;     // (J-1) * phi_L

        // store integrand in evaluations array
        evaluationsArrayPressure[samplingPointIndex][dofIndex] = integrand;
      }  // L

    }  // function evaluations

    // integrate all values for result vector entries at once
    EvaluationsDisplacementsType integratedValuesDisplacements = QuadratureDD::computeIntegral(evaluationsArrayDisplacements);
    EvaluationsPressureType integratedValuesPressure = QuadratureDD::computeIntegral(evaluationsArrayPressure);

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


    // add entries in result vector for pressure
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

  // assemble result vector
  if (communicateGhosts)
  {
    combinedVecResidual_->finishGhostManipulation();
  }
  //combinedVecResidual_->startGhostManipulation();
  //combinedVecResidual_->zeroGhostBuffer();
  //combinedVecResidual_->finishGhostManipulation();

  // now, solverVariableResidual_, which is the globalValues() of combinedVecResidual_, contains δW_int
}

template<typename Term,typename MeshType, int nDisplacementComponents>
void HyperelasticitySolver<Term,MeshType,nDisplacementComponents>::
materialComputeResidual(double loadFactor)
{
  // This computes the residual, i.e. the nonlinear function to be solved.
  // Compute Wint - Wext in variable solverVariableResidual_.
  //  output is solverVariableResidual_, a normal Vec, no Dirichlet BC dofs, also accessible by combinedVecResidual_
  //  input is solverVariableSolution_, a normal Vec, the same values have already been assigned to this->data_.displacements() and this->data_.pressure() (!)
  // before this method, values of u, v and p get stored to the data object by setUVP(solverVariableSolution_);

  LOG(DEBUG) << "materialComputeResidual";

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

  materialComputeInternalVirtualWork(false);    // compute without communicating ghost values, because startGhostManipulation has been called

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
    ierr = VecAXPY(solverVariableResidual_, -loadFactor, externalVirtualWorkDead_); CHKERRV(ierr);

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
    ierr = VecAXPY(solverVariableResidual_, -loadFactor, externalVirtualWorkDead_); CHKERRV(ierr);

    if (outputFiles)
    {
      static int evaluationNo = 0;  // counter how often this function was called

      // dumpVector(std::string filename, std::string format, Vec &vector, MPI_Comm mpiCommunicator, int componentNo=0, int nComponents=1);
      std::stringstream filename;
      filename << "out/id" << std::setw(3) << std::setfill('0') << evaluationNo++ << "r" << DihuContext::nRanksCommWorld();

      combinedVecResidual_->dumpGlobalNatural(filename.str());
    }
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
}

template<typename Term,typename MeshType, int nDisplacementComponents>
void HyperelasticitySolver<Term,MeshType,nDisplacementComponents>::
materialComputeExternalVirtualWorkDead()
{
  // compute δW_ext,dead = int_Ω B^L * phi^L * phi^M * δu^M dx + int_∂Ω T^L * phi^L * phi^M * δu^M dS

  LOG(DEBUG) << "materialComputeExternalVirtualWorkDead";

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
void HyperelasticitySolver<Term,MeshType,nDisplacementComponents>::
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

  typedef std::array<double, nDisplacementsDofsPerElement*nDisplacementsDofsPerElement> EvaluationsUVType;   // for vectorized types, phi_L*phi_M is the same for all elements, therefore double and not double_v_t
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

    // initialize top right and bottom left sub matrices to zero
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
  //ierr = MatDiagonalSet(combinedMatrixJacobian_->valuesGlobal(), zeros_, INSERT_VALUES); CHKERRV(ierr);

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

      Tensor2_v_t<D> deformationGradient = this->computeDeformationGradient(displacementsValues, inverseJacobianMaterial, xi);    // F
      double_v_t deformationGradientDeterminant;    // J
      Tensor2_v_t<D> inverseDeformationGradient = MathUtility::computeInverse(deformationGradient, deformationGradientDeterminant);  // F^-1

      Tensor2_v_t<D> rightCauchyGreen = this->computeRightCauchyGreenTensor(deformationGradient);  // C = F^T*F

      double_v_t rightCauchyGreenDeterminant;   // J^2
      Tensor2_v_t<D> inverseRightCauchyGreen = MathUtility::computeSymmetricInverse(rightCauchyGreen, rightCauchyGreenDeterminant);  // C^-1

      // fiber direction
      Vec3_v_t fiberDirection = displacementsFunctionSpace->template interpolateValueInElement<3>(elementalDirectionValues, xi);

      // invariants
      std::array<double_v_t,5> invariants = this->computeInvariants(rightCauchyGreen, rightCauchyGreenDeterminant, fiberDirection);  // I_1, I_2, I_3
      std::array<double_v_t,5> reducedInvariants = this->computeReducedInvariants(invariants, deformationGradientDeterminant); // Ibar_1, Ibar_2

      // pressure is the separately interpolated pressure for mixed formulation
      double_v_t pressure = pressureFunctionSpace->interpolateValueInElement(pressureValuesCurrentElement, xi);

      // Pk2 stress tensor S = S_vol + S_iso (p.234)
      //! compute 2nd Piola-Kirchhoff stress tensor S = 2*dPsi/dC and the fictitious PK2 Stress Sbar
      Tensor2_v_t<D> fictitiousPK2Stress;   // Sbar
      Tensor2_v_t<D> pk2StressIsochoric;    // S_iso
      Tensor2_v_t<D> pK2Stress = this->computePK2Stress(pressure, rightCauchyGreen, inverseRightCauchyGreen, reducedInvariants, deformationGradientDeterminant, fiberDirection,
                                                      fictitiousPK2Stress, pk2StressIsochoric
                                                    );

      std::array<Vec3,nDisplacementsDofsPerElement> gradPhi = displacementsFunctionSpace->getGradPhi(xi);
      // (column-major storage) gradPhi[L][a] = dphi_L / dxi_a
      // gradPhi[column][row] = gradPhi[dofIndex][i] = dphi_dofIndex/dxi_i, columnIdx = dofIndex, rowIdx = which direction


      Tensor4_v_t<D> elasticityTensor;
      Tensor4_v_t<D> fictitiousElasticityTensor;
      Tensor4_v_t<3> elasticityTensorIso;
      computeElasticityTensor(rightCauchyGreen, inverseRightCauchyGreen, deformationGradientDeterminant, pressure, reducedInvariants, fictitiousPK2Stress, pk2StressIsochoric, fiberDirection,
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
        LOG(FATAL) << "Deformation gradient " << deformationGradient << " has zero or negative determinant " << deformationGradientDeterminant
          << std::endl << "Geometry values in element " << elementNoLocal << ": " << geometryReferenceValues << std::endl
          << "Displacements at xi " << xi << ": " << displacementsValues;
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
              evaluationsArrayDisplacements[samplingPointIndex][index] = integrand;

            }  // b, bComponent
          }   // M, bDof
        }  // a, aComponent
      }  // L, aDof

      // add contributions of submatrix up and pu (lower left and upper right)

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
            evaluationsArrayPressure[samplingPointIndex][index] = integrand;

          }  // a
        }  // M
      }  // L

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
            evaluationsArrayUV[samplingPointIndex][index] = integrand;
          }   // M, mDof
        }   // L, lDof

      }  // if dynamic problem
    }   // sampling points

    // integrate all values for result vector entries at once
    EvaluationsDisplacementsType integratedValuesDisplacements = QuadratureDD::computeIntegral(evaluationsArrayDisplacements);
    EvaluationsPressureType integratedValuesPressure = QuadratureDD::computeIntegral(evaluationsArrayPressure);
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

    // add entries in result stiffness matrix for pressure (lower left and upper right parts, symmetric)
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
          const double integratedValue = integratedValuesUV[index];

          // integratedValue is only ∫_Ω ρ0 ϕ^L ϕ^M dV,
          // but we need 1/dt δ_ab ∫_Ω ρ0 ϕ^L ϕ^M dV

          for (int aComponent = 0; aComponent < D; aComponent++)           // a
          {
            for (int bComponent = 0; bComponent < D; bComponent++)           // b
            {
              if (aComponent != bComponent)
                continue;

              double resultingValue = 1./this->timeStepWidth_ * integratedValue;

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
}

template<typename Term,typename MeshType, int nDisplacementComponents>
template<typename double_v_t>
Tensor2<3,double_v_t> HyperelasticitySolver<Term,MeshType,nDisplacementComponents>::
computeDeformationGradient(const std::array<VecD<3,double_v_t>,DisplacementsFunctionSpace::nDofsPerElement()> &displacements,
                           const Tensor2<3,double_v_t> &inverseJacobianMaterial,
                           const std::array<double, 3> xi
                          )
{
  // compute the deformation gradient x_i,j = δ_ij + u_i,j
  // where j is dimensionColumn and i is component of the used Vec3's


  VLOG(3) << "compute deformation gradient Fij, displacements: " << displacements;

  const int D = 3;
  Tensor2<D,double_v_t> deformationGradient;

  const int nDofsPerElement = DisplacementsFunctionSpace::nDofsPerElement();

  // loop over dimension, i.e. columns of deformation gradient, j
  for (int dimensionColumn = 0; dimensionColumn < 3; dimensionColumn++)
  {
    // compute du_i/dX_dimensionColumn for all i at once
    VecD<3,double_v_t> du_dX({0});    // vector contains entries for all i, du_dXj with j = dimensionColumn

    VLOG(3) << " j = " << dimensionColumn;

    // loop over dof no., M
    for (int dofIndex = 0; dofIndex < nDofsPerElement; dofIndex++)
    {
      VLOG(3) << "  M = " << dofIndex;

      // compute dphi_dofIndex/dX_dimensionColumn
      double_v_t dphi_dX = 0;
      for (int l = 0; l < 3; l++)
      {
        VLOG(3) << "   l = " << l;
        double_v_t dphi_dxil = DisplacementsFunctionSpace::dphi_dxi(dofIndex, l, xi);
        double_v_t dxil_dX = inverseJacobianMaterial[dimensionColumn][l];     // inverseJacobianMaterial[j][l] = J_lj = dxi_l/dX_j

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

template<typename Term,typename MeshType, int nDisplacementComponents>
template<typename double_v_t>
Tensor2<3,double_v_t> HyperelasticitySolver<Term,MeshType,nDisplacementComponents>::
computeDeformationGradientTimeDerivative(const std::array<VecD<3,double_v_t>,DisplacementsFunctionSpace::nDofsPerElement()> &velocities,
                                         const Tensor2<3,double_v_t> &inverseJacobianMaterial,
                                         const std::array<double, 3> xi
                                        )
{
  // compute the time derivative of the deformation gradient d(x_i,j)/dt = d/dt(δ_ij + u_i,j) = v_i,j
  // where j is dimensionColumn and i is component of the used Vec3's


  VLOG(3) << "compute deformation gradient time derivative Fdot_ij, velocities: " << velocities;

  const int D = 3;
  Tensor2<D,double_v_t> deformationGradientTimeDerivative;

  const int nDofsPerElement = DisplacementsFunctionSpace::nDofsPerElement();

  // loop over dimension, i.e. columns of deformation gradient, j
  for (int dimensionColumn = 0; dimensionColumn < 3; dimensionColumn++)
  {
    // compute du_i/dX_dimensionColumn for all i at once
    VecD<3,double_v_t> du_dX({0});    // vector contains entries for all i, du_dXj with j = dimensionColumn

    VLOG(3) << " j = " << dimensionColumn;

    // loop over dof no., M
    for (int dofIndex = 0; dofIndex < nDofsPerElement; dofIndex++)
    {
      VLOG(3) << "  M = " << dofIndex;

      // compute dphi_dofIndex/dX_dimensionColumn
      double_v_t dphi_dX = 0;
      for (int l = 0; l < 3; l++)
      {
        VLOG(3) << "   l = " << l;
        double_v_t dphi_dxil = DisplacementsFunctionSpace::dphi_dxi(dofIndex, l, xi);
        double_v_t dxil_dX = inverseJacobianMaterial[dimensionColumn][l];     // inverseJacobianMaterial[j][l] = J_lj = dxi_l/dX_j

        VLOG(3) << "     dphi_dxil = " << dphi_dxil << ", dxil_dX = " << dxil_dX;

        // multiply dphi/dxi with dxi/dX to obtain dphi/dX
        dphi_dX += dphi_dxil * dxil_dX;
      }
      VLOG(3) << "   dphi_dX = " << dphi_dX;

      VLOG(3) << "   displ_M: " << velocities[dofIndex];
      du_dX += dphi_dX * velocities[dofIndex];   // vector-valued addition
    }
    VLOG(3) << " du_dXj: " << du_dX;

    deformationGradientTimeDerivative[dimensionColumn] = du_dX;
  }

  VLOG(3) << "deformationGradientTimeDerivative: " << deformationGradientTimeDerivative;
  return deformationGradientTimeDerivative;
}

template<typename Term,typename MeshType, int nDisplacementComponents>
template<typename double_v_t>
Tensor2<3,double_v_t> HyperelasticitySolver<Term,MeshType,nDisplacementComponents>::
computeRightCauchyGreenTensor(const Tensor2<3,double_v_t> &deformationGradient)
{
  // compute C = F^T*F where F is the deformationGradient and C is the right Cauchy-Green Tensor
  // the quantities are 3x3 tensors for the 3D case
  Tensor2<3,double_v_t> rightCauchyGreenTensor({std::array<double_v_t,3>({0})});

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

template<typename Term,typename MeshType, int nDisplacementComponents>
template<typename double_v_t>
std::array<double_v_t,5> HyperelasticitySolver<Term,MeshType,nDisplacementComponents>::
computeInvariants(const Tensor2<3,double_v_t> &rightCauchyGreen, const double_v_t rightCauchyGreenDeterminant, const VecD<3,double_v_t> fiberDirection)
{
  std::array<double_v_t,5> invariants;

  // I1 = tr(C)
  invariants[0] = 0.0;
  for (int i = 0; i < 3; i++)
  {
    invariants[0] += rightCauchyGreen[i][i];
  }

  // tr(C^2)
  double_v_t traceC2 = 0;
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

  // for a transversely isotropic material that also has 4th and 5th invariant
  if (Term::usesFiberDirection)
  {
    // I4 = a0 • C a0;
    double_v_t a0Ca0 = 0;
    for (int i = 0; i < 3; i++)
    {
      double_v_t ca0_i = 0;
      for (int j = 0; j < 3; j++)
      {
        ca0_i += rightCauchyGreen[j][i] * fiberDirection[j];
      }
      a0Ca0 += fiberDirection[i] * ca0_i;
    }

    invariants[3] = a0Ca0;

    // I5 = a0 • C^2 a0;

    double_v_t a0C2a0 = 0;
    for (int i = 0; i < 3; i++)
    {
      double_v_t c2a0_i = 0;
      for (int j = 0; j < 3; j++)
      {
        // compute C^2
        double_v_t c2_ij = 0;
        for (int k = 0; k < 3; k++)
        {
          // C^2_ij = C_ik * C_kj
          c2_ij += rightCauchyGreen[k][i] * rightCauchyGreen[j][k];
        }
        c2a0_i += c2_ij * fiberDirection[j];
      }
      a0C2a0 += fiberDirection[i] * c2a0_i;
    }

    invariants[4] = a0C2a0;
    //LOG(DEBUG) << "computed I4: " << invariants[3] << ", I5: " << invariants[4] << ", a0: " << fiberDirection;
  }

  return invariants;
}

template<typename Term,typename MeshType, int nDisplacementComponents>
template<typename double_v_t>
std::array<double_v_t,5> HyperelasticitySolver<Term,MeshType,nDisplacementComponents>::
computeReducedInvariants(const std::array<double_v_t,5> invariants, const double_v_t deformationGradientDeterminant)
{
  std::array<double_v_t,5> reducedInvariants;

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

  reducedInvariants[0] = MathUtility::pow(deformationGradientDeterminant, factor23) * invariants[0];
  reducedInvariants[1] = MathUtility::pow(deformationGradientDeterminant, factor43) * invariants[1];

  // for a transversely isotropic material that also has 4th and 5th invariant
  if (Term::usesFiberDirection)
  {
    reducedInvariants[2] = 1;  // not used, because I3 = det C = 1 constant
    reducedInvariants[3] = MathUtility::pow(deformationGradientDeterminant, factor23) * invariants[3];
    reducedInvariants[4] = MathUtility::pow(deformationGradientDeterminant, factor43) * invariants[4];
  }

  return reducedInvariants;
}

template<typename Term,typename MeshType, int nDisplacementComponents>
template<typename double_v_t>
Tensor2<3,double_v_t> HyperelasticitySolver<Term,MeshType,nDisplacementComponents>::
computePK2Stress(const double_v_t pressure,                             //< [in] pressure value p
                 const Tensor2<3,double_v_t> &rightCauchyGreen,         //< [in] C
                 const Tensor2<3,double_v_t> &inverseRightCauchyGreen,  //< [in] C^{-1}
                 const std::array<double_v_t,5> reducedInvariants,      //< [in] the reduced invariants Ibar_1, Ibar_2
                 const double_v_t deformationGradientDeterminant,       //< [in] J = det(F)
                 VecD<3,double_v_t> fiberDirection,                     //< [in] a0, direction of fibers
                 Tensor2<3,double_v_t> &fictitiousPK2Stress,            //< [out] Sbar, the fictitious 2nd Piola-Kirchhoff stress tensor
                 Tensor2<3,double_v_t> &pk2StressIsochoric              //< [out] S_iso, the isochoric part of the 2nd Piola-Kirchhoff stress tensor
                )
{
  // compute the PK2 stress tensor as S=2*dPsi/dC
  // for explanation see pdf document
  auto dPsi_dIbar1Expression = SEMT::deriv_t(Term::strainEnergyDensityFunctionIsochoric, Term::Ibar1);
  auto dPsi_dIbar2Expression = SEMT::deriv_t(Term::strainEnergyDensityFunctionIsochoric, Term::Ibar2);

  std::vector<double_v_t> reducedInvariantsVector(reducedInvariants.begin(), reducedInvariants.end());

  const double_v_t dPsi_dIbar1 = ExpressionHelper<double_v_t>::apply(dPsi_dIbar1Expression, reducedInvariantsVector);
  const double_v_t dPsi_dIbar2 = ExpressionHelper<double_v_t>::apply(dPsi_dIbar2Expression, reducedInvariantsVector);

  const double_v_t Ibar1 = reducedInvariants[0];
  const double_v_t J = deformationGradientDeterminant;

  double factor23 = -2./3;
  double_v_t factor1 = 2*(dPsi_dIbar1 + Ibar1*dPsi_dIbar2);
  double_v_t factor2 = -2*dPsi_dIbar2;
  double_v_t factorJ23 = MathUtility::pow(J, factor23);

  Tensor2<3,double_v_t> pK2Stress;

  // compute fictitiousPK2Stress:
  // Sbar = factor1*I + factor2*Cbar, Cbar = J^{-2/3}*C

  // row index
  for (int i=0; i<3; i++)
  {
    // column index
    for (int j=0; j<3; j++)
    {
      int delta_ij = (i == j? 1 : 0);
      const double_v_t cBar = factorJ23 * rightCauchyGreen[j][i];
      fictitiousPK2Stress[j][i] = factor1 * delta_ij + factor2 * cBar;

      //if (i == j)
      //LOG(DEBUG) << "  C_" << i << j << " = " << rightCauchyGreen[j][i] << ", J=" << J << ", factorJ23 = " << factorJ23 << ", SBar_" << i << j << ": " << fictitiousPK2Stress[j][i]
      //     << " = " << factor1 << " * " << delta_ij << " + " << factor2 << "*" << cBar;
    }
  }

  //LOG(DEBUG) << "in computePK2Stress, factor1: " << factor1 << ", factor2: " << factor2 << ", C: " << rightCauchyGreen << ", Sbar: " << fictitiousPK2Stress << ", reducedInvariants: " << reducedInvariants;

  // add term for 4th invariant
  if (Term::usesFiberDirection)
  {
    auto dPsi_dIbar4Expression = SEMT::deriv_t(Term::strainEnergyDensityFunctionIsochoric, Term::Ibar4);
    const double_v_t dPsi_dIbar4 = ExpressionHelper<double_v_t>::apply(dPsi_dIbar4Expression, reducedInvariantsVector);
    double_v_t factor4 = 2*dPsi_dIbar4;

    for (int i=0; i<3; i++)     // alternative indices: A
    {
      // column index
      for (int j=0; j<3; j++)     // alternative indices: B
      {
        fictitiousPK2Stress[j][i] += factor4 * fiberDirection[i] * fiberDirection[j];
      }
    }

    //LOG(DEBUG) << "fiberDirection: " << fiberDirection << ", factor4: " << factor4
    //  << ", dPsi_dIbar4Expression: " << dPsi_dIbar4Expression << ", invariants: " << reducedInvariantsVector;
  }

  // add term for 5th invariant (not implemented)

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
      const double_v_t sVol = J * pressure * inverseRightCauchyGreen[j][i];         // S_vol = J * p * C^{-1}_AB

      // compute P : Sbar
      double_v_t pSbar = 0;
      // row index
      for (int k=0; k<3; k++)        // alternative indices: C
      {
        const int delta_ik = (i == k? 1 : 0);

        // column index
        for (int l=0; l<3; l++)            // alternative indices: D
        {
          const int delta_jl = (j == l? 1 : 0);

          // this is a non-symmetric version for Ii but it is also correct, the symmetric version would be given by Ii = (δ_AC*δ_BD + δ_AD*δ_BC) / 2
          const int Ii = delta_ik * delta_jl;
          const double_v_t Cc = inverseRightCauchyGreen[j][i] * rightCauchyGreen[l][k];     // CC = C^{-1}_AB * C_CD
          const double_v_t Pp = (Ii - 1./3 * Cc);

          //LOG(DEBUG) << "    PP_" << i << j << k << l << " = " << Pp << " = " << Ii << "-1/3*" << Cc << " (" << inverseRightCauchyGreen[j][i] << "," << rightCauchyGreen[l][k] << "), Ii: " << Ii;

          pSbar += Pp * fictitiousPK2Stress[l][k];
        }
      }


      // isochoric stress
      const double_v_t sIso = factorJ23 * pSbar;
      pk2StressIsochoric[j][i] = sIso;
      //LOG(DEBUG) << "    PSbar_" << i << j << ": " << pSbar << ", J^-2/3: " << factorJ23 << ", Siso_" << i << j << ": " << sIso;

      // total stress is sum of volumetric and isochoric part, sVol = J*p*C^{-1}_AB, sIso = j^{-2/3}*(Ii-1/3*Cc)*Sbar
      pK2Stress[j][i] = sVol + sIso;

      // debugging
      //pK2Stress[j][i] = sIso;

      //VLOG(2) << "set pk2Stress_" << i << j << " = " << pK2Stress[j][i];

      //if (i == j)
      //  LOG(DEBUG) << "  ccs: " << ccs << " C:Sbar: " << cSbar << ", factorJ23: " << factorJ23 << ", Svol_" << i << j << " = " << sVol << ", Siso_" << i << j << " = " << sIso << ", S = " << pK2Stress[j][i];
    }  // j
  }  // i

  //for debugging, check symmetry of PK2 stress and if it is correct according to Mooney-Rivlin formula
  if (VLOG_IS_ON(2))
  {
    const double_v_t errorTolerance = 1e-14;
    bool pK2IsSymmetric = true;
    for (int a=0; a<3; a++)
    {
      for (int b=0; b<3; b++)
      {
        if (Vc::any_of(MathUtility::abs(pK2Stress[a][b] - pK2Stress[b][a]) > errorTolerance))
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
    const double_v_t factorJ23 = MathUtility::pow(deformationGradientDeterminant, factor23);

    const double c0 = SEMT::Parameter<0>::get_value();
    const double c1 = SEMT::Parameter<1>::get_value();
    //const double c0 = PARAM(0).get_value();   //< material parameter
    //const double c1 = PARAM(1).get_value();   //< material parameter

    const double_v_t Ibar1 = reducedInvariants[0];

    const double_v_t gamma1 = 2*(c0 + c1*Ibar1);
    const double_v_t gamma2 = -2*c1;

    bool mismatch = false;
    for (int a = 0; a < 3; a++)
    {
      for (int b = 0; b < 3; b++)
      {
        const int delta_ab = (a == b? 1 : 0);
        double_v_t cBar = factorJ23 * rightCauchyGreen[b][a];
        double_v_t sBar = gamma1*delta_ab + gamma2*cBar;

        if (Vc::any_of(MathUtility::abs(sBar - fictitiousPK2Stress[b][a]) > errorTolerance))
        {
          LOG(ERROR) << "mismatch in Sbar_" << a << b << ": derived: " << fictitiousPK2Stress << ", Mooney Rivlin explicit formula: " << sBar;
          mismatch = true;
        }
      }
    }
    if (!mismatch)
      LOG_N_TIMES(2,DEBUG) << "Sbar is correct!";
  }

  return pK2Stress;
}

template<typename Term,typename MeshType, int nDisplacementComponents>
void HyperelasticitySolver<Term,MeshType,nDisplacementComponents>::
computePK2StressField()
{
  //LOG(TRACE) << "computePK2StressField";

  //this->data_.pK2Stress()->startGhostManipulation();
  this->data_.pK2Stress()->zeroGhostBuffer();
  this->data_.deformationGradient()->zeroGhostBuffer();
  this->data_.deformationGradientTimeDerivative()->zeroGhostBuffer();

  // get pointer to function space
  std::shared_ptr<DisplacementsFunctionSpace> displacementsFunctionSpace = this->data_.displacementsFunctionSpace();
  std::shared_ptr<PressureFunctionSpace> pressureFunctionSpace = this->data_.pressureFunctionSpace();

  const int D = 3;  // dimension
  const int nDisplacementsDofsPerElement = DisplacementsFunctionSpace::nDofsPerElement();
  const int nPressureDofsPerElement = PressureFunctionSpace::nDofsPerElement();
  const int nElementsLocal = displacementsFunctionSpace->nElementsLocal();

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

    // get displacements field values for element
    std::array<Vec3_v_t,nDisplacementsDofsPerElement> velocitiesValues;
    this->data_.velocities()->getElementValues(elementNoLocalv, velocitiesValues);

    std::array<double_v_t,nPressureDofsPerElement> pressureValuesCurrentElement;
    this->data_.pressure()->getElementValues(elementNoLocalv, pressureValuesCurrentElement);

    std::array<Vec3_v_t,nDisplacementsDofsPerElement> elementalDirectionValues;
    this->data_.fiberDirection()->getElementValues(elementNoLocalv, elementalDirectionValues);

    // get indices of element-local dofs
    std::array<dof_no_v_t,27> dofNosLocal = this->displacementsFunctionSpace_->getElementDofNosLocal(elementNoLocalv);

    //LOG(DEBUG) << "el " << elementNoLocal << ", geometryRef: " << geometryReferenceValues << ", displacements: " << displacementsValues << ", p: " << pressureValuesCurrentElement;

    // loop over nodes of this element
    for (int elementalNodeNo = 0; elementalNodeNo < 27; elementalNodeNo++)
    {
      dof_no_v_t dofNoLocal = dofNosLocal[elementalNodeNo];

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
      Tensor2_v_t<D> jacobianMaterial = DisplacementsFunctionSpace::computeJacobian(geometryReferenceValues, xi);
      double_v_t jacobianDeterminant;
      Tensor2_v_t<D> inverseJacobianMaterial = MathUtility::computeInverse(jacobianMaterial, jacobianDeterminant);

      // jacobianMaterial[columnIdx][rowIdx] = dX_rowIdx/dxi_columnIdx
      // inverseJacobianMaterial[columnIdx][rowIdx] = dxi_rowIdx/dX_columnIdx because of inverse function theorem

      // F
      Tensor2_v_t<D> deformationGradient = this->computeDeformationGradient(displacementsValues, inverseJacobianMaterial, xi);

      double_v_t deformationGradientDeterminant = MathUtility::computeDeterminant(deformationGradient);  // J

      // compute Fdot values
      Tensor2_v_t<D> Fdot = computeDeformationGradientTimeDerivative(velocitiesValues, inverseJacobianMaterial, xi);

      // store F values
      std::array<double_v_t,9> deformationGradientValues;
      std::array<double_v_t,9> deformationGradientTimeDerivativeValues;
      for (int j = 0; j < 3; j++)
      {
        for (int i = 0; i < 3; i++)
        {
          // row-major
          deformationGradientValues[j*3+i] = deformationGradient[j][i];
          deformationGradientTimeDerivativeValues[j*3+i] = Fdot[j][i];
        }
      }
      this->data_.deformationGradient()->setValue(dofNoLocal, deformationGradientValues, INSERT_VALUES);
      this->data_.deformationGradientTimeDerivative()->setValue(dofNoLocal, deformationGradientTimeDerivativeValues, INSERT_VALUES);


      Tensor2_v_t<D> rightCauchyGreen = this->computeRightCauchyGreenTensor(deformationGradient);  // C = F^T*F

      double_v_t rightCauchyGreenDeterminant;   // J^2
      Tensor2_v_t<D> inverseRightCauchyGreen = MathUtility::computeSymmetricInverse(rightCauchyGreen, rightCauchyGreenDeterminant);  // C^-1

      // fiber direction
      Vec3_v_t fiberDirection = displacementsFunctionSpace->template interpolateValueInElement<3>(elementalDirectionValues, xi);

      // invariants
      std::array<double_v_t,5> invariants = this->computeInvariants(rightCauchyGreen, rightCauchyGreenDeterminant, fiberDirection);  // I_1, I_2, I_3, I_4, I_5
      std::array<double_v_t,5> reducedInvariants = this->computeReducedInvariants(invariants, deformationGradientDeterminant); // Ibar_1, Ibar_2, Ibar_4, Ibar_5

      // pressure is the separately interpolated pressure for mixed formulation
      double_v_t pressure = pressureFunctionSpace->interpolateValueInElement(pressureValuesCurrentElement, xi);

      // Pk2 stress tensor S = S_vol + S_iso (p.234)
      //! compute 2nd Piola-Kirchhoff stress tensor S = 2*dPsi/dC and the fictitious PK2 Stress Sbar
      Tensor2_v_t<D> fictitiousPK2Stress;   // Sbar
      Tensor2_v_t<D> pk2StressIsochoric;    // S_iso
      Tensor2_v_t<D> pK2Stress = this->computePK2Stress(pressure, rightCauchyGreen, inverseRightCauchyGreen, reducedInvariants, deformationGradientDeterminant, fiberDirection,
                                                      fictitiousPK2Stress, pk2StressIsochoric
                                                    );


      std::array<double_v_t,6> valuesInVoigtNotation({pK2Stress[0][0], pK2Stress[1][1], pK2Stress[2][2], pK2Stress[0][1], pK2Stress[1][2], pK2Stress[0][2]});

      //LOG(DEBUG) << "node " << dofNoLocal << " pk2: " << valuesInVoigtNotation;
      this->data_.pK2Stress()->setValue(dofNoLocal, valuesInVoigtNotation, INSERT_VALUES);
    }
  }
  this->data_.pK2Stress()->zeroGhostBuffer();
  this->data_.pK2Stress()->finishGhostManipulation();
  this->data_.pK2Stress()->startGhostManipulation();

  this->data_.deformationGradient()->zeroGhostBuffer();
  this->data_.deformationGradient()->finishGhostManipulation();
  this->data_.deformationGradient()->startGhostManipulation();

  this->data_.deformationGradientTimeDerivative()->zeroGhostBuffer();
  this->data_.deformationGradientTimeDerivative()->finishGhostManipulation();
  this->data_.deformationGradientTimeDerivative()->startGhostManipulation();
}

//! compute the material elasticity tensor
template<typename Term,typename MeshType, int nDisplacementComponents>
template<typename double_v_t>
void HyperelasticitySolver<Term,MeshType,nDisplacementComponents>::
computeElasticityTensor(const Tensor2<3,double_v_t> &rightCauchyGreen,         //< [in] C
                        const Tensor2<3,double_v_t> &inverseRightCauchyGreen,  //< [in] C^{-1}
                        double_v_t deformationGradientDeterminant,             //< [in] J = det(F)
                        double_v_t pressure,                                   //< [in] pressure value p
                        std::array<double_v_t,5> reducedInvariants,            //< [in] the reduced invariants Ibar_1, Ibar_2
                        const Tensor2<3,double_v_t> &fictitiousPK2Stress,      //< [in] Sbar, the fictitious 2nd Piola-Kirchhoff stress tensor
                        const Tensor2<3,double_v_t> &pk2StressIsochoric,       //< [in] S_iso, the isochoric part of the 2nd Piola-Kirchhoff stress tensor
                        VecD<3,double_v_t> fiberDirection,                     //< [in] a0, direction of fibers
                        Tensor4<3,double_v_t> &fictitiousElasticityTensor,     //< [out] fictitious Elasticity tensor CCbar_{ABCD}
                        Tensor4<3,double_v_t> &elasticityTensorIso,            //< [out] CCiso_{ABCD}
                        Tensor4<3,double_v_t> &elasticityTensor                //< [out] elasticity tensor CC_{ABCD}
                       )
{
  // compute the elasticity tensor as CC=2*dS(C)/dC
  // for explanation see pdf document
  const int D = 3;

  // compute preliminary variables that are independent of the indices a,b,c,d
  auto dPsi_dIbar1Expression = SEMT::deriv_t(Term::strainEnergyDensityFunctionIsochoric, Term::Ibar1);
  auto dPsi_dIbar2Expression = SEMT::deriv_t(Term::strainEnergyDensityFunctionIsochoric, Term::Ibar2);

  auto d2Psi_dIbar1Ibar1Expression = SEMT::deriv_t(dPsi_dIbar1Expression, Term::Ibar1);
  auto d2Psi_dIbar1Ibar2Expression = SEMT::deriv_t(dPsi_dIbar1Expression, Term::Ibar2);
  auto d2Psi_dIbar2Ibar2Expression = SEMT::deriv_t(dPsi_dIbar2Expression, Term::Ibar2);

  std::vector<double_v_t> reducedInvariantsVector(reducedInvariants.begin(), reducedInvariants.end());

  //const double dPsi_dIbar1 = dPsi_dIbar1Expression.apply(reducedInvariantsVector);
  const double_v_t dPsi_dIbar2       = ExpressionHelper<double_v_t>::apply(dPsi_dIbar2Expression, reducedInvariantsVector);
  const double_v_t d2Psi_dIbar1Ibar1 = ExpressionHelper<double_v_t>::apply(d2Psi_dIbar1Ibar1Expression, reducedInvariantsVector);
  const double_v_t d2Psi_dIbar1Ibar2 = ExpressionHelper<double_v_t>::apply(d2Psi_dIbar1Ibar2Expression, reducedInvariantsVector);
  const double_v_t d2Psi_dIbar2Ibar2 = ExpressionHelper<double_v_t>::apply(d2Psi_dIbar2Ibar2Expression, reducedInvariantsVector);

  const double_v_t Ibar1 = reducedInvariants[0];

  // formula for C_iso: Holzapfel "Nonlinear Solid Mechanics" p.255
  // formula for C_bar: Holzapfel "Nonlinear Solid Mechanics" p.262
  // compute factors for Cbar
  const double_v_t factor1 = 4*(d2Psi_dIbar1Ibar1 + 2*Ibar1*d2Psi_dIbar1Ibar2 + dPsi_dIbar2 + MathUtility::sqr(Ibar1)*d2Psi_dIbar2Ibar2);
  const double_v_t factor2 = -4*(d2Psi_dIbar1Ibar2 + Ibar1*d2Psi_dIbar2Ibar2);
  const double_v_t factor3 = 4*d2Psi_dIbar2Ibar2;
  const double_v_t factor4 = -4*dPsi_dIbar2;

  if (false)
  {
    LOG(DEBUG) << "elasticity tensor, Ψ: " << Term::strainEnergyDensityFunctionIsochoric;
    LOG(DEBUG) << "∂Ψ/∂Ibar1: " << dPsi_dIbar1Expression;
    LOG(DEBUG) << "∂Ψ/∂Ibar2: " << dPsi_dIbar2Expression << " = " << dPsi_dIbar2;
    LOG(DEBUG) << "∂2Ψ/(∂Ibar1 ∂Ibar1): " << d2Psi_dIbar1Ibar1Expression << " = " << d2Psi_dIbar1Ibar1;
    LOG(DEBUG) << "∂2Ψ/(∂Ibar1 ∂Ibar1): " << d2Psi_dIbar1Ibar2Expression << " = " << d2Psi_dIbar1Ibar2;
    LOG(DEBUG) << "∂2Ψ/(∂Ibar1 ∂Ibar1): " << d2Psi_dIbar2Ibar2Expression << " = " << d2Psi_dIbar2Ibar2;
    LOG(DEBUG) << "factor1: " << factor1;
    LOG(DEBUG) << "factor2: " << factor2;
    LOG(DEBUG) << "factor3: " << factor3;
    LOG(DEBUG) << "factor4: " << factor4;
  }

  const double_v_t J = deformationGradientDeterminant;

  const double_v_t factorJ23 = MathUtility::pow(J,-2./3);   // J^{-2/3}
  const double_v_t factorJ43 = MathUtility::pow(J,-4./3);   // J^{-4/3}

  // distinct entries, only those have to be computed as the rest is symmetric
  const int indices[21][4] = {
    {0,0,0,0},{0,1,0,0},{0,2,0,0},{1,1,0,0},{1,2,0,0},{2,2,0,0},{0,1,0,1},{0,2,0,1},{1,1,0,1},{1,2,0,1},
    {2,2,0,1},{0,2,0,2},{1,1,0,2},{1,2,0,2},{2,2,0,2},{1,1,1,1},{1,2,1,1},{2,2,1,1},{1,2,1,2},{2,2,1,2},
    {2,2,2,2}
  };

  // to compute all entries and verify the symmetry, use the following code (needs further adjustments at  the loop boundaries: "entryNo<21" -> "entryNo<81")
#if 0
  int indices[81][4];

  int entryNo = 0;
  for (int a = 0; a < D; a++)
  {
    for (int b = 0; b < D; b++)
    {
      for (int c = 0; c < D; c++)
      {
        for (int d = 0; d < D; d++, entryNo++)
        {
          indices[entryNo][0] = a;
          indices[entryNo][1] = b;
          indices[entryNo][2] = c;
          indices[entryNo][3] = d;
        }
      }
    }
  }
  assert(entryNo == 81);
#endif

  // loop over distinct entries in elasticity tensor
  for (int entryNo = 0; entryNo < 21; entryNo++)
  {
    // get indices of current entry
    const int a = indices[entryNo][0];
    const int b = indices[entryNo][1];
    const int c = indices[entryNo][2];
    const int d = indices[entryNo][3];

    const double_v_t cInv_ab = inverseRightCauchyGreen[b][a];
    const double_v_t cInv_cd = inverseRightCauchyGreen[d][c];

    // compute C_vol
    const double_v_t cInvDotCInv = 0.5 * (inverseRightCauchyGreen[c][a] * inverseRightCauchyGreen[d][b] + inverseRightCauchyGreen[d][a] * inverseRightCauchyGreen[c][b]);

    const double_v_t cInvDyadCInv = inverseRightCauchyGreen[b][a] * inverseRightCauchyGreen[d][c];
    const double_v_t jp = deformationGradientDeterminant * pressure;

    const double_v_t cVol = jp * cInvDyadCInv - 2 * jp * cInvDotCInv;

    // compute C_iso
    //                          ab  ef     gh   cd
    // compute contribution from  P : Cbar : P^T
    double_v_t pCbarPT_abcd = 0.0;

    for (int g = 0; g < D; g++)
    {
      for (int h = 0; h < D; h++)
      {
        int delta_gh = (g == h? 1 : 0);

        double_v_t pCbar_abgh = 0;
        const double_v_t c_gh = rightCauchyGreen[h][g];

        for (int e = 0; e < D; e++)
        {
          int delta_eh = (e == h? 1 : 0);
          int delta_eg = (e == g? 1 : 0);

          for (int f = 0; f < D; f++)
          {
            int delta_fg = (f == g? 1 : 0);
            int delta_fh = (f == h? 1 : 0);
            int delta_ef = (e == f? 1 : 0);
            const double_v_t c_ef = rightCauchyGreen[f][e];

            int iI_efgh = delta_eg * delta_fh;
            int iIbar_efgh = delta_eh * delta_fg;

            // symmetric version
            double_v_t sS_efgh = 0.5 * (iI_efgh + iIbar_efgh);

            const double_v_t summand1 = factor1 * delta_ef * delta_gh;
            const double_v_t summand2 = factor2 * (delta_ef * factorJ23*c_gh + factorJ23*c_ef * delta_gh);
            const double_v_t summand3 = factor3 * factorJ23*c_ef * factorJ23*c_gh;
            const double_v_t summand4 = factor4 * sS_efgh;

            double_v_t sum = summand1 + summand2 + summand3 + summand4;

            // terms for 4th invariant
            if (Term::usesFiberDirection)
            {
              auto dPsi_dIbar4Expression = SEMT::deriv_t(Term::strainEnergyDensityFunctionIsochoric, Term::Ibar4);

              auto d2Psi_dIbar1Ibar4Expression = SEMT::deriv_t(dPsi_dIbar1Expression, Term::Ibar4);
              auto d2Psi_dIbar2Ibar4Expression = SEMT::deriv_t(dPsi_dIbar2Expression, Term::Ibar4);
              auto d2Psi_dIbar4Ibar4Expression = SEMT::deriv_t(dPsi_dIbar4Expression, Term::Ibar4);

              const double_v_t d2Psi_dIbar1Ibar4 = ExpressionHelper<double_v_t>::apply(d2Psi_dIbar1Ibar4Expression, reducedInvariantsVector);
              const double_v_t d2Psi_dIbar2Ibar4 = ExpressionHelper<double_v_t>::apply(d2Psi_dIbar2Ibar4Expression, reducedInvariantsVector);
              const double_v_t d2Psi_dIbar4Ibar4 = ExpressionHelper<double_v_t>::apply(d2Psi_dIbar4Ibar4Expression, reducedInvariantsVector);

              const double_v_t factor5 = 4*(d2Psi_dIbar1Ibar4 + Ibar1 * d2Psi_dIbar2Ibar4);
              const double_v_t factor6 = -4*d2Psi_dIbar2Ibar4;
              const double_v_t factor7 = 4*d2Psi_dIbar4Ibar4;

              const double_v_t summand5 = factor5 * (delta_ef       * fiberDirection[g] * fiberDirection[h] + fiberDirection[e] * fiberDirection[f] * delta_gh);
              const double_v_t summand6 = factor6 * (factorJ23*c_ef * fiberDirection[g] * fiberDirection[h] + fiberDirection[e] * fiberDirection[f] * factorJ23*c_gh);
              const double_v_t summand7 = factor7 * fiberDirection[e] * fiberDirection[f] * fiberDirection[g] * fiberDirection[h];

              sum += summand5 + summand6 + summand7;

              // terms for 5th invariant are not implemented, this means that the strain energy function can only depend on Ibar4
            }

            const double_v_t ccBar_efgh = factorJ43 * sum;

            //LOG(DEBUG) << "     CCBar_" << e << f << g << h << ": " << factorJ43 << "*(" << summand1 << "+" << summand2 << "+" << summand3 << "+" << summand4 << "), summand4 = " << factor4 << "*" << sS_efgh;

            fictitiousElasticityTensor[h][g][f][e] = ccBar_efgh;

            int delta_ae = (a == e? 1 : 0);
            int delta_bf = (b == f? 1 : 0);

            const double_v_t iI_abef = delta_ae * delta_bf;
            const double_v_t p_abef = iI_abef - 1./D*cInv_ab*c_ef;

            pCbar_abgh += p_abef * ccBar_efgh;

          }  // f
        }  // e

        int delta_cg = (c == g? 1 : 0);
        int delta_dh = (d == h? 1 : 0);
        const double_v_t iI_cdgh = delta_cg * delta_dh;

        const double_v_t p_cdgh = iI_cdgh - 1./D*cInv_cd*c_gh;     // (P^T)_ghcd = P_cdgh

        const double_v_t pT_ghcd = p_cdgh;

        pCbarPT_abcd += pCbar_abgh * pT_ghcd;

      }  // h
    }  // g

    // compute contribution from  2/3*J^{-2/3}*Sbar : C P_tilde
    double_v_t sBarC_abcd = 0.;
    for (int g = 0; g < D; g++)
    {
      for (int h = 0; h < D; h++)
      {
        sBarC_abcd += fictitiousPK2Stress[h][g] * rightCauchyGreen[h][g];
      }  // h
    }  // g

    const double_v_t pTilde_abcd = cInvDotCInv - 1./3 * cInvDyadCInv;

    const double_v_t sBarCP_abcd = 2./3. * factorJ23 * sBarC_abcd * pTilde_abcd;

    // compute contribution from -2./3*(CInv dyad Siso + Siso dyad CInv)
    const double_v_t cInvSiso = -2./3 * (cInv_ab * pk2StressIsochoric[d][c] + pk2StressIsochoric[b][a] * cInv_cd);

    // compute C_iso
    const double_v_t cIso = pCbarPT_abcd + sBarCP_abcd + cInvSiso;

    //LOG(DEBUG) << a<<b<<c<<d<<": compute CCIso: " << cIso << ": " << pCbarPT_abcd << "," << sBarCP_abcd << "," << cInvSiso;

#if 0
    // store entry in CCIso, this is only returned from this method for debugging purposes, only for use in materialTesting
    // therefore it can be omitted here
    elasticityTensorIso[d][c][b][a] = cIso;
    elasticityTensorIso[c][d][b][a] = cIso;
    elasticityTensorIso[d][c][a][b] = cIso;
    elasticityTensorIso[c][d][a][b] = cIso;
    elasticityTensorIso[b][a][d][c] = cIso;
    elasticityTensorIso[a][b][d][c] = cIso;
    elasticityTensorIso[b][a][c][d] = cIso;
    elasticityTensorIso[a][b][c][d] = cIso;
#endif
    // compute C
    const double_v_t c_abcd = cVol + cIso;

    // debugging
    //const double_v_t c_abcd = cIso;

    // store entry, store all symmetric entries at once, because tensor has major and minor symmetries (CC_abcd = CC_cdab and CC_abcd = CC_bacd)
    elasticityTensor[d][c][b][a] = c_abcd;
    elasticityTensor[c][d][b][a] = c_abcd;
    elasticityTensor[d][c][a][b] = c_abcd;
    elasticityTensor[c][d][a][b] = c_abcd;
    elasticityTensor[b][a][d][c] = c_abcd;
    elasticityTensor[a][b][d][c] = c_abcd;
    elasticityTensor[b][a][c][d] = c_abcd;
    elasticityTensor[a][b][c][d] = c_abcd;
  }

  // verify major and minor symmetries of elasticity tensor
  if (VLOG_IS_ON(1))
  {
  //LOG_N_TIMES(2,DEBUG) << "elasticityTensor: " << elasticityTensor;

    const double_v_t errorTolerance = 1e-12;
    std::string name("C");
    for (int a = 0; a < 3; a++)
    {
      for (int b = 0; b < 3; b++)
      {
        for (int c = 0; c < 3; c++)
        {
          for (int d = 0; d < 3; d++)
          {
            if (Vc::any_of(MathUtility::abs(elasticityTensor[a][b][c][d] - elasticityTensor[b][a][c][d]) > errorTolerance))
            {
              LOG(ERROR) << name << "[" <<a<< "][" <<b<< "][" << c<< "][" << d<< "] != " << name << "[" <<b<< "][" <<a<< "][" << c<< "][" << d<< "] (" <<elasticityTensor[b][a][c][d]<< " != " <<elasticityTensor[a][b][c][d]<< ") - minor symmetry violated" << std::endl;
            }
            if (Vc::any_of(MathUtility::abs(elasticityTensor[a][b][c][d] - elasticityTensor[a][b][d][c]) > errorTolerance))
            {
              LOG(ERROR) << name << "[" <<a<< "][" <<b<< "][" << c<< "][" << d<< "] != " << name << "[" <<a<< "][" <<b<< "][" << d<< "][" << c<< "] (" <<elasticityTensor[a][b][d][c]<< " != " <<elasticityTensor[a][b][c][d]<< ") - minor symmetry violated" << std::endl;
            }
            if (Vc::any_of(MathUtility::abs(elasticityTensor[a][b][c][d] - elasticityTensor[c][d][a][b]) > errorTolerance))
            {
              LOG(ERROR) << name << "[" <<a<< "][" <<b<< "][" << c<< "][" << d<< "] != " << name << "[" << c<< "][" << d<< "][" <<a<< "][" <<b<< "] (" <<elasticityTensor[c][d][a][b]<< " != " <<elasticityTensor[a][b][c][d]<< ") - major symmetry violated" << std::endl;
            }
          }
        }
      }
    }
    LOG_N_TIMES(2,DEBUG) << "elasticity tensor checked for minor symmetries";
    //LOG(DEBUG) << "elasticity tensor checked for minor symmetries";
  }
}

template<typename Term,typename MeshType, int nDisplacementComponents>
template<typename double_v_t>
Tensor2<3,double_v_t> HyperelasticitySolver<Term,MeshType,nDisplacementComponents>::
computePSbar(const Tensor2<3,double_v_t> &fictitiousPK2Stress, const Tensor2<3,double_v_t> &rightCauchyGreen)
{
  // only needed for debugging in materialTesting
  double_v_t determinant;
  Tensor2<3,double_v_t> inverseRightCauchyGreen = MathUtility::computeInverse<3>(rightCauchyGreen, determinant);

  //ab cd
  //  P : Sbar
  Tensor2<3,double_v_t> PSbar;
  for (int a = 0; a < 3; a++)
  {
    for (int b = 0; b < 3; b++)
    {
      double_v_t cInv_ab = inverseRightCauchyGreen[b][a];

      double_v_t psbar_ab = 0;
      for (int c = 0; c < 3; c++)
      {
        int delta_ac = (a == c? 1 : 0);
        for (int d = 0; d < 3; d++)
        {
          int delta_bd = (b == d? 1 : 0);
          const double_v_t ii_abcd = delta_ac * delta_bd;

          double_v_t c_cd = rightCauchyGreen[d][c];

          double_v_t p_abcd = ii_abcd - 1./3*cInv_ab * c_cd;

          //LOG(DEBUG) << "  .PP_" << a << b << c << d << ": " << p_abcd << "=" << ii_abcd << "-1/3*" << cInv_ab * c_cd << " (" << cInv_ab << "," << c_cd << "), Ii: " << ii_abcd;

          const double_v_t sBar_cd = fictitiousPK2Stress[d][c];
          psbar_ab += p_abcd * sBar_cd;
        }
      }


      //LOG(DEBUG) << " PSBar_" << a << b << ": " << psbar_ab;
      PSbar[b][a] = psbar_ab;
    }
  }
  return PSbar;
}

} // namespace
