#include "spatial_discretization/finite_element_method/solid_mechanics/solid_mechanics_common.h"

#include <Python.h>  // has to be the first included header
#include <memory>
#include <array>
#include <petscksp.h>
#include <petscsys.h>
#include <cmath>

#include "semt/Semt.h"
#include "semt/Shortcuts.h"
#include "control/types.h"
#include "spatial_discretization/finite_element_method/solid_mechanics/elasticity_tensor.h"
#include "basis_on_mesh/00_basis_on_mesh_base_dim.h"
#include "utility/python_utility.h"
#include "utility/matrix.h"

#ifdef QUADRATURE_TEST
#include "quadrature/gauss.h"
#include "quadrature/clenshaw_curtis.h"
#include "quadrature/newton_cotes.h"
#endif

namespace SpatialDiscretization
{


template<typename BasisOnMeshType,typename BasisOnMeshTypeForUtility, typename QuadratureType, typename Term>
void SolidMechanicsCommon<BasisOnMeshType,BasisOnMeshTypeForUtility,QuadratureType,Term>::
setStiffnessMatrix()
{
  // call setStiffnessMatrix of the derived class. That, in turn, calls setStiffnessMatrixEntriesForDisplacements of this class.
  this->setStiffnessMatrix(PETSC_NULL);
}

template<typename BasisOnMeshType,typename BasisOnMeshTypeForUtility, typename MixedQuadratureType, typename Term>
void SolidMechanicsCommon<BasisOnMeshType,BasisOnMeshTypeForUtility,MixedQuadratureType,Term>::
setStiffnessMatrixEntriesForDisplacements(Mat tangentStiffnessMatrix)
{
  typedef typename BasisOnMeshType::HighOrderBasisOnMesh BasisOnMesh;

  // get pointer to mesh object
  std::shared_ptr<BasisOnMesh> mesh = this->data_.mesh();

  const int D = BasisOnMesh::dim();  // = 2 or 3
  const int nDofsPerElement = BasisOnMesh::nDofsPerElement();
  const int nElements = mesh->nElements();
  const int nUnknowsPerElement = nDofsPerElement*D;    // 3 directions for displacements per dof

  // define shortcuts for quadrature
  typedef typename MixedQuadratureType::HighOrderQuadrature QuadratureType;
  typedef Quadrature::TensorProduct<D,QuadratureType> QuadratureDD;

  // define type to hold evaluations of integrand for stiffness matrix
  typedef MathUtility::Matrix<nUnknowsPerElement,nUnknowsPerElement> EvaluationsType;
  typedef std::array<
            EvaluationsType,
            QuadratureDD::numberEvaluations()
          > EvaluationsArrayType;     // evaluations[nGP^D](nUnknows,nUnknows)

  // setup arrays used for integration
  std::array<VecD<D>, QuadratureDD::numberEvaluations()> samplingPoints = QuadratureDD::samplingPoints();
  EvaluationsArrayType evaluationsArray;

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

  //std::array<std::array<double,nUnknowsPerElement*nUnknowsPerElement>,QuadratureExactDD::numberEvaluations()> evaluationsExactArray;
  EvaluationsExactArrayType evaluationsExactArray;
#endif

  PetscErrorCode ierr;

  LOG(TRACE) << "setStiffnessMatrix - compute tangent stiffness matrix";
  LOG(DEBUG) << "nUnknowsPerElement: " << nUnknowsPerElement<<", n evaluations for quadrature: " << QuadratureDD::numberEvaluations();

  // initialize values to zero
  if (!tangentStiffnessMatrixInitialized_)
  {
    // there are no previous non-zero entries, so loop over all later needed entries and manually set to zero

    int cntr = 1;
    // loop over elements
    for (element_no_t elementNo = 0; elementNo < mesh->nElements(); elementNo++)
    {
      auto dofNo = mesh->getElementDofNos(elementNo);

      for (int aDof = 0; aDof < nDofsPerElement; aDof++)
      {
        for (int aComponent = 0; aComponent < D; aComponent++)
        {
          for (int bDof = 0; bDof < nDofsPerElement; bDof++)
          {
            for (int bComponent = 0; bComponent < D; bComponent++)
            {
              dof_no_t matrixRowIndex = dofNo[aDof]*D + aComponent;
              dof_no_t matrixColumnIndex = dofNo[bDof]*D + bComponent;

              VLOG(3) << " initialize tangentStiffnessMatrix (("<<aDof<<","<<aComponent<<"),("<<bDof<<","<<bComponent<<")), dofs (" << dofNo[aDof] << ","<<dofNo[bDof]<<"), entry ( " << matrixRowIndex << "," << matrixColumnIndex << ") (no. " << cntr++ << ")";
              ierr = MatSetValue(tangentStiffnessMatrix, matrixRowIndex, matrixColumnIndex, 0.0, INSERT_VALUES); CHKERRV(ierr);
            }
          }
        }
      }
    }
  }
  else
  {
    // zero all previous non-zero entries, keeps matrix structure
    MatZeroEntries(tangentStiffnessMatrix);

    VLOG(3) << " zero tangentStiffnessMatrix";
  }

  // loop over elements
  for (int elementNo = 0; elementNo < nElements; elementNo++)
  {

#ifdef QUADRATURE_TEST
    Control::PerformanceMeasurement::start("stiffnessMatrixDisplacements");
#endif

    // get geometry field of reference configuration, the geometry field is always 3D, also in 2D problems
    std::array<Vec3,nDofsPerElement> geometryReferenceValues;
    this->data_.geometryReference().getElementValues(elementNo, geometryReferenceValues);

    // get displacement field
    std::array<VecD<D>,nDofsPerElement> displacementValues;
    this->data_.displacements().getElementValues(elementNo, displacementValues);

    // For mixed formulation get the pressure values of this element. This is done only in the derived class for mixed formulation and not in the derived class for penalty formulation.
    this->preparePressureInterpolation(elementNo);

    // loop over integration points (e.g. gauss points) for displacement field
    for (unsigned int samplingPointIndex = 0; samplingPointIndex < samplingPoints.size(); samplingPointIndex++)
    {
      // get parameter values of current sampling point
      VecD<D> xi = samplingPoints[samplingPointIndex];

      // compute the 3xD jacobian of the parameter space to world space mapping
      Tensor2<D> jacobianMaterial = MathUtility::transformToDxD<D,D>(BasisOnMesh::computeJacobian(geometryReferenceValues, xi));
      double jacobianDeterminant;
      Tensor2<D> inverseJacobianMaterial = MathUtility::computeInverse<D>(jacobianMaterial, jacobianDeterminant);
      // jacobianMaterial[columnIdx][rowIdx] = dX_rowIdx/dxi_columnIdx
      // inverseJacobianMaterial[columnIdx][rowIdx] = dxi_rowIdx/dX_columnIdx because of inverse function theorem

      // F
      Tensor2<D> deformationGradient = this->computeDeformationGradient(displacementValues, inverseJacobianMaterial, xi);
      double deformationGradientDeterminant = MathUtility::computeDeterminant<D>(deformationGradient);  // J

      Tensor2<D> rightCauchyGreen = this->computeRightCauchyGreenTensor(deformationGradient);  // C = F^T*F

      double rightCauchyGreenDeterminant;   // J^2
      Tensor2<D> inverseRightCauchyGreen = MathUtility::computeSymmetricInverse<D>(rightCauchyGreen, rightCauchyGreenDeterminant);  // C^-1

      // invariants
      std::array<double,3> invariants = this->computeInvariants(rightCauchyGreen, rightCauchyGreenDeterminant);  // I_1, I_2, I_3
      std::array<double,2> reducedInvariants = this->computeReducedInvariants(invariants, deformationGradientDeterminant); // Ibar_1, Ibar_2

      // pressure is the artificialPressure for penalty formulation or the separately interpolated pressure for mixed formulation
      double artificialPressureTilde;
      const double artificialPressure = this->getPressure(deformationGradientDeterminant, xi, artificialPressureTilde);

      // Pk2 stress tensor S = S_vol + S_iso (p.234)
      Tensor2<D> fictitiousPK2Stress;
      Tensor2<D> pk2StressIsochoric;
      Tensor2<D> PK2Stress = this->computePK2Stress(artificialPressure, rightCauchyGreen, inverseRightCauchyGreen, reducedInvariants, deformationGradientDeterminant, fictitiousPK2Stress, pk2StressIsochoric);
      // elasticity tensor C_{ijkl}
      ElasticityTensor elasticity = this->computeElasticityTensor(artificialPressure, artificialPressureTilde, rightCauchyGreen, inverseRightCauchyGreen, fictitiousPK2Stress, pk2StressIsochoric, deformationGradientDeterminant, reducedInvariants);

      std::array<VecD<D>,nDofsPerElement> gradPhi = mesh->getGradPhi(xi);
      // (column-major storage) gradPhi[M][a] = dphi_M / dxi_a
      // gradPhi[column][row] = gradPhi[dofIndex][i] = dphi_dofIndex/dxi_i, columnIdx = dofIndex, rowIdx = which direction

      VLOG(2) << "geometryReferenceValues: " << geometryReferenceValues;
      VLOG(2) << "Jacobian: J_phi=" << jacobianMaterial;
      VLOG(2) << "jacobianDeterminant: J=" << jacobianDeterminant;
      VLOG(2) << "inverseJacobianMaterial: J_phi^-1=" << inverseJacobianMaterial;
      VLOG(2) << "deformationGradient: F=" << deformationGradient;
      VLOG(2) << "deformationGradientDeterminant: det F=" << deformationGradientDeterminant;
      VLOG(2) << "rightCauchyGreen: C=" << rightCauchyGreen;
      VLOG(2) << "rightCauchyGreenDeterminant: det C=" << rightCauchyGreenDeterminant;
      VLOG(2) << "inverseRightCauchyGreen: C^-1=" << inverseRightCauchyGreen;
      VLOG(2) << "invariants: I1,I2,I3: " << invariants;
      VLOG(2) << "reducedInvariants: Ibar1, Ibar2: " << reducedInvariants;
      VLOG(2) << "artificialPressure: p=" << artificialPressure << ", artificialPressureTilde: pTilde=" << artificialPressureTilde;
      VLOG(2) << "fictitiousPK2Stress: Sbar=" << fictitiousPK2Stress;
      VLOG(2) << "pk2StressIsochoric: S_iso=" << pk2StressIsochoric;
      VLOG(2) << "PK2Stress: S=" << PK2Stress;
      VLOG(2) << "gradPhi: " << gradPhi;
      VLOG(2) << "elasticity: C=" << elasticity;

      if (deformationGradientDeterminant < 1e-12)   // if any entry of the deformation gradient is negative
      {
        LOG(FATAL) << "Deformation gradient " << deformationGradient << " has zero or negative determinant " << deformationGradientDeterminant;
      }

      // loop over pairs of basis functions and evaluate integrand at xi
      for (int aDof = 0; aDof < nDofsPerElement; aDof++)           // index over dofs, each dof has D components, L in derivation
      {
        for (int aComponent = 0; aComponent < D; aComponent++)     // lower-case a in derivation, index over displacement components
        {
          // compute index of degree of freedom and component (matrix row index)
          const int i = D*aDof + aComponent;

          for (int bDof = 0; bDof < nDofsPerElement; bDof++)       // index M in derivation
          {
            for (int bComponent = 0; bComponent < D; bComponent++)     // lower-case b in derivation, index over displacement components
            {

              // compute index of degree of freedom and component (index of entry in vector of unknows / matrix column index)
              const int j = D*bDof + bComponent;

              // compute integrand[i][j]
              double integrand = 0.0;
              for (int dInternal = 0; dInternal < D; dInternal++)     // capital D in derivation
              {
                for (int bInternal = 0; bInternal < D; bInternal++)     // capital B in derivation
                {
                  // ----------------------------
                  // compute derivatives of phi
                  double dphiL_dXB = 0.0;
                  double dphiM_dXD = 0.0;

                  // helper index k for multiplication with inverse Jacobian
                  for (int k = 0; k < D; k++)
                  {
                    // compute dphiL/dXB from dphiL/dxik and dxik/dXB
                    const double dphiL_dxik = gradPhi[aDof][k];    // dphi_L/dxik
                    const double dxik_dXB = inverseJacobianMaterial[bInternal][k];  // inverseJacobianMaterial[B][k] = J_kB = dxi_k/dX_B

                    dphiL_dXB += dphiL_dxik * dxik_dXB;

                    // compute dphiM/dXD from dphiM/dxik and dxik/dXD
                    const double dphiM_dxik = gradPhi[bDof][k];    // dphi_M/dxik
                    const double dxik_dXD = inverseJacobianMaterial[dInternal][k];  // inverseJacobianMaterial[D][k] = J_kD = dxi_k/dX_D

                    dphiM_dXD += dphiM_dxik * dxik_dXD;
                  }

                  // ----------------------------
                  // compute ktilde_abBD
                  const int delta_ab = (aComponent == bComponent? 1 : 0);
                  const double S_BD = PK2Stress[dInternal][bInternal];

                  // compute ffC = FaA * FbC * C_ABCD
                  double ffC = 0.0;
                  for (int aInternal = 0; aInternal < D; aInternal++)     // capital A in derivation
                  {
                    for (int cInternal = 0; cInternal < D; cInternal++)    // capital C in derivation
                    {
                      const double FaA = deformationGradient[aInternal][aComponent];
                      const double FbC = deformationGradient[cInternal][bComponent];
                      const double cc = elasticity.getEntry(aInternal, bInternal, cInternal, dInternal);
                      ffC += FaA * FbC * cc;
                    }
                  }

                  const double ktilde_abBD = delta_ab*S_BD + ffC;

                  // k_ij = int dphi_La/dX_B * dphi_Mb/dX_D * ktilde_abBD
                  integrand += dphiL_dXB * ktilde_abBD * dphiM_dXD;
                }
              }

              VLOG(2) << "    (L,a),(M,b)=(" << aDof << ","<<aComponent<<"),("<<bDof<<","<<bComponent<<"), (i,j)=("<<i<<","<<j<<"), integrand=" << integrand;

              // store integrand in evaluations array
              evaluationsArray[samplingPointIndex](i,j) = integrand * fabs(jacobianDeterminant);
            }  // b
          }  // M
        }  // a
      }  // L

    }  // function evaluations

    // integrate all values for the (i,j) dof pairs at once
    EvaluationsType integratedValues = QuadratureDD::computeIntegral(evaluationsArray);

#ifdef QUADRATURE_TEST
    Control::PerformanceMeasurement::stop("stiffnessMatrixDisplacements");
#endif

#ifdef QUADRATURE_TEST
    // loop over integration points (e.g. gauss points) for displacement field
    for (unsigned int samplingPointIndex = 0; samplingPointIndex < samplingPointsExact.size(); samplingPointIndex++)
    {
      // get parameter values of current sampling point
      VecD<D> xi = samplingPointsExact[samplingPointIndex];

      // compute the 3xD jacobian of the parameter space to world space mapping
      Tensor2<D> jacobianMaterial = MathUtility::transformToDxD<D,D>(BasisOnMesh::computeJacobian(geometryReferenceValues, xi));
      double jacobianDeterminant;
      Tensor2<D> inverseJacobianMaterial = MathUtility::computeInverse<D>(jacobianMaterial, jacobianDeterminant);
      // jacobianMaterial[columnIdx][rowIdx] = dX_rowIdx/dxi_columnIdx
      // inverseJacobianMaterial[columnIdx][rowIdx] = dxi_rowIdx/dX_columnIdx because of inverse function theorem

      // F
      Tensor2<D> deformationGradient = this->computeDeformationGradient(displacementValues, inverseJacobianMaterial, xi);
      double deformationGradientDeterminant = MathUtility::computeDeterminant<D>(deformationGradient);  // J

      Tensor2<D> rightCauchyGreen = this->computeRightCauchyGreenTensor(deformationGradient);  // C = F^T*F

      double rightCauchyGreenDeterminant;   // J^2
      Tensor2<D> inverseRightCauchyGreen = MathUtility::computeSymmetricInverse<D>(rightCauchyGreen, rightCauchyGreenDeterminant);  // C^-1

      // invariants
      std::array<double,3> invariants = this->computeInvariants(rightCauchyGreen, rightCauchyGreenDeterminant);  // I_1, I_2, I_3
      std::array<double,2> reducedInvariants = this->computeReducedInvariants(invariants, deformationGradientDeterminant); // Ibar_1, Ibar_2

      // pressure is the artificialPressure for penalty formulation or the separately interpolated pressure for mixed formulation
      double artificialPressureTilde;
      const double artificialPressure = this->getPressure(deformationGradientDeterminant, xi, artificialPressureTilde);

      // Pk2 stress tensor S = S_vol + S_iso (p.234)
      Tensor2<D> fictitiousPK2Stress;
      Tensor2<D> pk2StressIsochoric;
      Tensor2<D> PK2Stress = this->computePK2Stress(artificialPressure, rightCauchyGreen, inverseRightCauchyGreen, reducedInvariants, deformationGradientDeterminant, fictitiousPK2Stress, pk2StressIsochoric);
      // elasticity tensor C_{ijkl}
      ElasticityTensor elasticity = this->computeElasticityTensor(artificialPressure, artificialPressureTilde, rightCauchyGreen, inverseRightCauchyGreen, fictitiousPK2Stress, pk2StressIsochoric, deformationGradientDeterminant, reducedInvariants);

      std::array<VecD<D>,nDofsPerElement> gradPhi = mesh->getGradPhi(xi);
      // (column-major storage) gradPhi[M][a] = dphi_M / dxi_a
      // gradPhi[column][row] = gradPhi[dofIndex][i] = dphi_dofIndex/dxi_i, columnIdx = dofIndex, rowIdx = which direction

      if (deformationGradientDeterminant < 1e-12)   // if any entry of the deformation gradient is negative
      {
        LOG(FATAL) << "Deformation gradient " << deformationGradient << " has zero or negative determinant " << deformationGradientDeterminant;
      }

      // loop over pairs of basis functions and evaluate integrand at xi
      for (int aDof = 0; aDof < nDofsPerElement; aDof++)           // index over dofs, each dof has D components, L in derivation
      {
        for (int aComponent = 0; aComponent < D; aComponent++)     // lower-case a in derivation, index over displacement components
        {
          // compute index of degree of freedom and component (matrix row index)
          const int i = D*aDof + aComponent;

          for (int bDof = 0; bDof < nDofsPerElement; bDof++)       // index M in derivation
          {
            for (int bComponent = 0; bComponent < D; bComponent++)     // lower-case b in derivation, index over displacement components
            {

              // compute index of degree of freedom and component (index of entry in vector of unknows / matrix column index)
              const int j = D*bDof + bComponent;

              // compute integrand[i][j]
              double integrand = 0.0;
              for (int dInternal = 0; dInternal < D; dInternal++)     // capital D in derivation
              {
                for (int bInternal = 0; bInternal < D; bInternal++)     // capital B in derivation
                {
                  // ----------------------------
                  // compute derivatives of phi
                  double dphiL_dXB = 0.0;
                  double dphiM_dXD = 0.0;

                  // helper index k for multiplication with inverse Jacobian
                  for (int k = 0; k < D; k++)
                  {
                    // compute dphiL/dXB from dphiL/dxik and dxik/dXB
                    const double dphiL_dxik = gradPhi[aDof][k];    // dphi_L/dxik
                    const double dxik_dXB = inverseJacobianMaterial[bInternal][k];  // inverseJacobianMaterial[B][k] = J_kB = dxi_k/dX_B

                    dphiL_dXB += dphiL_dxik * dxik_dXB;

                    // compute dphiM/dXD from dphiM/dxik and dxik/dXD
                    const double dphiM_dxik = gradPhi[bDof][k];    // dphi_M/dxik
                    const double dxik_dXD = inverseJacobianMaterial[dInternal][k];  // inverseJacobianMaterial[D][k] = J_kD = dxi_k/dX_D

                    dphiM_dXD += dphiM_dxik * dxik_dXD;
                  }

                  // ----------------------------
                  // compute ktilde_abBD
                  const int delta_ab = (aComponent == bComponent? 1 : 0);
                  const double S_BD = PK2Stress[dInternal][bInternal];

                  // compute ffC = FaA * FbC * C_ABCD
                  double ffC = 0.0;
                  for (int aInternal = 0; aInternal < D; aInternal++)     // capital A in derivation
                  {
                    for (int cInternal = 0; cInternal < D; cInternal++)    // capital C in derivation
                    {
                      const double FaA = deformationGradient[aInternal][aComponent];
                      const double FbC = deformationGradient[cInternal][bComponent];
                      const double cc = elasticity.getEntry(aInternal, bInternal, cInternal, dInternal);
                      ffC += FaA * FbC * cc;
                    }
                  }

                  const double ktilde_abBD = delta_ab*S_BD + ffC;

                  // k_ij = int dphi_La/dX_B * dphi_Mb/dX_D * ktilde_abBD
                  integrand += dphiL_dXB * ktilde_abBD * dphiM_dXD;
                }
              }

              //LOG(DEBUG) << "(L,a),(M,b)=(" << aDof << ","<<aComponent<<"),("<<bDof<<","<<bComponent<<"), (i,j)=("<<i<<","<<j<<"), integrand=" << integrand;

              // store integrand in evaluations array
              evaluationsExactArray[samplingPointIndex](i,j) = integrand * fabs(jacobianDeterminant);
//              evaluationsExactArray[samplingPointIndex][i*nUnknowsPerElement + j] = integrand * fabs(jacobianDeterminant);

            }  // b
          }  // M
        }  // a
      }  // L

    }  // function evaluations

    // integrate all values for result vector entries at once
    EvaluationsType integratedValuesExact = QuadratureExactDD::computeIntegral(evaluationsExactArray);

    Control::PerformanceMeasurement::measureError("stiffnessMatrixDisplacements", (integratedValues - integratedValuesExact)/integratedValuesExact);
#endif

    // get indices of element-local dofs
    std::array<dof_no_t,nDofsPerElement> dofNo = mesh->getElementDofNos(elementNo);

    // set entries in stiffness matrix
    // loop over indices of unknows ((aDof,aComponent), (bDof,bComponent)), i.e. (i,j)
    for (int aDof = 0; aDof < nDofsPerElement; aDof++)
    {
      for (int aComponent = 0; aComponent < D; aComponent++)
      {
        // compute index of degree of freedom and component (matrix row index)
        const int i = D*aDof + aComponent;

        for (int bDof = 0; bDof < nDofsPerElement; bDof++)
        {
          for (int bComponent = 0; bComponent < D; bComponent++)
          {
            // compute index of degree of freedom and component (index of entry in vector of unknows / matrix column index)
            const int j = D*bDof + bComponent;

            // integrate value and set entry in stiffness matrix
            double integratedValue = integratedValues(i,j);

            // Compute indices in stiffness matrix. For each dof there are D values for the D displacement components.
            // Therefore the nDofsPerElement number is not the number of unknows.
            dof_no_t matrixRowIndex = dofNo[aDof]*D + aComponent;
            dof_no_t matrixColumnIndex = dofNo[bDof]*D + bComponent;

            //VLOG(2) << "  pair (("<<aDof<<","<<aComponent<<"),("<<bDof<<","<<bComponent<<")) = (" <<i<<","<<j<<"), dofs (" << dofNo[aDof] << ","<<dofNo[bDof]<<")";
            //VLOG(2) << "      matrix indices ("<<matrixRowIndex<<","<<matrixColumnIndex<<"), integrated value: "<<integratedValue;

            ierr = MatSetValue(tangentStiffnessMatrix, matrixRowIndex, matrixColumnIndex, integratedValue, ADD_VALUES); CHKERRV(ierr);

          }  // bComponent
        }  // bDof
      }  // aComponent
    }  // aDof
  }  // elementNo

}

template<typename BasisOnMeshType,typename BasisOnMeshTypeForUtility, typename QuadratureType, typename Term>
Mat &SolidMechanicsCommon<BasisOnMeshType,BasisOnMeshTypeForUtility,QuadratureType,Term>::
tangentStiffnessMatrix()
{
  return this->data_.tangentStiffnessMatrix();
}

template<typename BasisOnMeshType,typename BasisOnMeshTypeForUtility, typename QuadratureType, typename Term>
Vec &SolidMechanicsCommon<BasisOnMeshType,BasisOnMeshTypeForUtility,QuadratureType,Term>::
rightHandSide()
{
  return this->data_.rightHandSide().values();
}

template<typename BasisOnMeshType,typename BasisOnMeshTypeForUtility, typename QuadratureType, typename Term>
Vec &SolidMechanicsCommon<BasisOnMeshType,BasisOnMeshTypeForUtility,QuadratureType,Term>::
displacements()
{
  return this->data_.displacements().values();
}

template<typename BasisOnMeshType,typename BasisOnMeshTypeForUtility, typename MixedQuadratureType, typename Term>
void SolidMechanicsCommon<BasisOnMeshType,BasisOnMeshTypeForUtility,MixedQuadratureType,Term>::
computeExternalVirtualWork(Vec &resultVec)
{
  LOG(TRACE) << "computeExternalVirtualWork";

  // get pointer to mesh object
  typedef typename BasisOnMeshType::HighOrderBasisOnMesh BasisOnMesh;  // for mixed formulation get the high order BasisOnMesh
  typedef typename MixedQuadratureType::HighOrderQuadrature QuadratureType;  // for mixed formulation get the high order quadrature

  std::shared_ptr<BasisOnMesh> mesh = this->data_.mesh();

  const int D = BasisOnMesh::dim();  // = 2 or 3
  //const int nDofs = mesh->nDofs();
  const int nDofsPerElement = BasisOnMesh::nDofsPerElement();
  //const int nElements = mesh->nElements();
  // define shortcuts for quadrature
  typedef Quadrature::TensorProduct<D,QuadratureType> QuadratureDD;
  typedef Quadrature::TensorProduct<D-1,QuadratureType> QuadratureSurface;

  // define type to hold evaluations of integrand for result vector, volume integrals for body forces
  typedef std::array<double, nDofsPerElement*D> EvaluationsType;
  typedef std::array<
            EvaluationsType,
            QuadratureDD::numberEvaluations()
          > EvaluationsArrayType;     // evaluations[nGP^D](nDofs*D)

  // setup arrays used for integration
  std::array<VecD<D>, QuadratureDD::numberEvaluations()> samplingPoints = QuadratureDD::samplingPoints();
  EvaluationsArrayType evaluationsArray{};

  // define type to hold evaluations of integrand for result vector, surface integrals for traction
  typedef std::array<
            EvaluationsType,
            QuadratureSurface::numberEvaluations()
          > EvaluationsArraySurfaceType;     // evaluations[nGP^D](nDofs*D)

  std::array<std::array<double,D-1>, QuadratureSurface::numberEvaluations()> samplingPointsSurface = QuadratureSurface::samplingPoints();
  EvaluationsArraySurfaceType evaluationsArraySurface;

  PetscErrorCode ierr;

  // zero result vector
  ierr = VecZeroEntries(resultVec); CHKERRV(ierr);

  // -------------- body force in reference configuration --------------
  // loop over elements in bodyForceReferenceConfiguration_ that have a force vector specified
  for (typename std::vector<std::pair<element_no_t, VecD<D>>>::const_iterator iter = this->bodyForceReferenceConfiguration_.begin();
       iter != this->bodyForceReferenceConfiguration_.end(); iter++)
  {
    int elementNo = iter->first;
    VecD<D> forceVector = iter->second;

    // check if element no is valid
    if (elementNo < 0 || elementNo > mesh->nElements())
    {
      LOG(ERROR) << "Element " << elementNo << " for which body force (ref.conf.) is specified is invalid (number of elements: " << mesh->nElements() << ")";
      continue;
    }


    // loop over integration points (e.g. gauss points)
    for (unsigned int samplingPointIndex = 0; samplingPointIndex < samplingPoints.size(); samplingPointIndex++)
    {
      // get parameter values of current sampling point
      VecD<D> xi = samplingPoints[samplingPointIndex];

      // loop over dofs of element
      for (int dofIndex = 0; dofIndex < nDofsPerElement; dofIndex++)
      {
        VecD<D> dofIntegrand = forceVector * MathUtility::sqr(mesh->phi(dofIndex, xi));

        // store integrand in evaluations array
        for (int i = 0; i < D; i++)
        {
          evaluationsArray[samplingPointIndex][dofIndex*D+i] = dofIntegrand[i];
        }

      }  // dofIndex

    }  // function evaluations

    // integrate all values for result vector entries at once
    EvaluationsType integratedValues = QuadratureDD::computeIntegral(evaluationsArray);

    // get indices of element-local dofs
    std::array<dof_no_t,nDofsPerElement> dofNo = mesh->getElementDofNos(elementNo);

    // add entries in result vector
    // loop over indices of unknows (dofIndex,dofComponent)
    for (int dofIndex = 0; dofIndex < nDofsPerElement; dofIndex++)
    {
      for (int dofComponent = 0; dofComponent < D; dofComponent++)
      {
        // compute index of degree of freedom and component (integrade values vector index)
        const int i = D*dofIndex + dofComponent;

        // get integrated value
        double integratedValue = integratedValues[i];

        // Compute indices in result vector. For each dof there are D values for the D displacement components
        // Therefore the nDofsPerElement number is not the number of unknows.
        dof_no_t resultVectorIndex = dofNo[dofIndex]*D + dofComponent;

        ierr = VecSetValue(resultVec, resultVectorIndex, integratedValue, ADD_VALUES); CHKERRV(ierr);

      }  // dofComponent
    }  // dofIndex
  }  // elementNo

  // -------------- body force in current configuration --------------
  // loop over elements in bodyForceCurrentConfiguration_ that have a force vector specified
  for (typename std::vector<std::pair<element_no_t, VecD<D>>>::const_iterator iter = this->bodyForceCurrentConfiguration_.begin();
       iter != this->bodyForceCurrentConfiguration_.end(); iter++)
  {
    int elementNo = iter->first;
    VecD<D> forceVector = iter->second;

    // check if element no is valid
    if (elementNo < 0 || elementNo > mesh->nElements())
    {
      LOG(ERROR) << "Element " << elementNo << " for which body force (cur.conf.) is specified is invalid (number of elements: " << mesh->nElements() << ")";
      continue;
    }

    // get geometry field of reference configuration, this is always 3D
    std::array<Vec3,nDofsPerElement> geometryReferenceValues;
    this->data_.geometryReference().getElementValues(elementNo, geometryReferenceValues);

    // loop over integration points (e.g. gauss points)
    for (unsigned int samplingPointIndex = 0; samplingPointIndex < samplingPoints.size(); samplingPointIndex++)
    {
      // get parameter values of current sampling point
      VecD<D> xi = samplingPoints[samplingPointIndex];

      // compute the 3xD jacobian of the parameter space to world space mapping
      Tensor2<D> jacobianMaterial = MathUtility::transformToDxD<D,D>(BasisOnMesh::computeJacobian(geometryReferenceValues, xi));
      double jacobianDeterminant = MathUtility::computeDeterminant<D>(jacobianMaterial);

      // loop over dofs of element
      for (int dofIndex = 0; dofIndex < nDofsPerElement; dofIndex++)
      {
        // scale the force vector with the Jacobian determinant: B = J*b
        // for incompressible material this should have no effect, J=1

        VecD<D> dofIntegrand = forceVector * jacobianDeterminant * MathUtility::sqr(mesh->phi(dofIndex, xi));

        // store integrand in evaluations array
        for (int i = 0; i < D; i++)
        {
          evaluationsArray[samplingPointIndex][dofIndex*D + i] = dofIntegrand[i];
        }

      }  // dofIndex

    }  // function evaluations

    // integrate all values for result vector entries at once
    EvaluationsType integratedValues = QuadratureDD::computeIntegral(evaluationsArray);

    // get indices of element-local dofs
    std::array<dof_no_t,nDofsPerElement> dofNo = mesh->getElementDofNos(elementNo);

    // add entries in result vector
    // loop over indices of unknows (dofIndex,dofComponent)
    for (int dofIndex = 0; dofIndex < nDofsPerElement; dofIndex++)
    {
      for (int dofComponent = 0; dofComponent < D; dofComponent++)
      {
        // compute index of degree of freedom and component (integrade values vector index)
        const int i = D*dofIndex + dofComponent;

        // get integrated value
        double integratedValue = integratedValues[i];

        // Compute indices in result vector. For each dof there are D values for the D displacement components.
        // Therefore the nDofsPerElement number is not the number of unknows.
        dof_no_t resultVectorIndex = dofNo[dofIndex]*D + dofComponent;

        ierr = VecSetValue(resultVec, resultVectorIndex, integratedValue, ADD_VALUES); CHKERRV(ierr);

      }  // dofComponent
    }  // dofIndex
  }  // elementNo

  // -------------- surface traction in reference configuration --------------
  for (typename std::vector<typename SolidMechanicsBoundaryConditions<BasisOnMeshType,Term>::TractionBoundaryCondition>::const_iterator iter = this->tractionReferenceConfiguration_.begin();
       iter != this->tractionReferenceConfiguration_.end(); iter++)
  {
    //  element_no_t elementGlobalNo;
    //
    //  Mesh::face_t face;
    //  std::vector<std::pair<dof_no_t, Vec3>> dofVectors;  //<element-local dof no, value>

    element_no_t elementNo = iter->elementGlobalNo;

    // check if element no is valid
    if (elementNo < 0 || elementNo > mesh->nElements())
    {
      LOG(ERROR) << "Element " << elementNo << " for which traction (ref.conf.) is specified is invalid (number of elements: " << mesh->nElements() << ")";
      continue;
    }

    LOG(DEBUG) << "tractionReferenceConfiguration_ on element " << elementNo;

    // loop over integration points (e.g. gauss points)
    for (unsigned int samplingPointIndex = 0; samplingPointIndex < samplingPointsSurface.size(); samplingPointIndex++)
    {
      // get parameter values of current sampling point
      std::array<double,D-1> xiSurface = samplingPointsSurface[samplingPointIndex];
      VecD<D> xi = Mesh::getXiOnFace(iter->face, xiSurface);

      // set all entries to 0
      evaluationsArraySurface[samplingPointIndex] = {0.0};

      // compute the traction at xi by summing contributions from elemental dofs
      VecD<D> tractionReferenceConfiguration({0.0});
      for (typename std::vector<std::pair<dof_no_t, VecD<D>>>::const_iterator dofVectorsIter = iter->dofVectors.begin();
           dofVectorsIter != iter->dofVectors.end();
           dofVectorsIter++)
      {
        int dofIndex = dofVectorsIter->first;
        VecD<D> forceVector = dofVectorsIter->second;

        tractionReferenceConfiguration += forceVector * mesh->phi(dofIndex, xi);
      }

      // loop over dofs of element with given traction, those have potentially non-zero virtual external work contribution
      for (typename std::vector<std::pair<dof_no_t, VecD<D>>>::const_iterator dofVectorsIter = iter->dofVectors.begin();
           dofVectorsIter != iter->dofVectors.end();
           dofVectorsIter++)
      {
        int dofIndex = dofVectorsIter->first;
        VecD<D> forceVector = dofVectorsIter->second;

        VecD<D> dofIntegrand = tractionReferenceConfiguration * mesh->phi(dofIndex, xi);

        LOG(DEBUG) << "  dofIndex " << dofIndex << ", xi=" << xi << ", traction: " << tractionReferenceConfiguration[0]
          << " phi = " << mesh->phi(dofIndex, xi);

        // store integrand in evaluations array
        for (int i = 0; i < D; i++)
        {
          evaluationsArraySurface[samplingPointIndex][dofIndex*D+i] = dofIntegrand[i];
        }
      }  // dofIndex

    }  // function evaluations

    // integrate all values for result vector entries at once
    EvaluationsType integratedValues = QuadratureSurface::computeIntegral(evaluationsArraySurface);

    // get indices of element-local dofs
    std::array<dof_no_t,nDofsPerElement> dofNo = mesh->getElementDofNos(elementNo);

    // add entries in result vector
    // loop over indices of unknows (dofIndex,dofComponent)
    for (int dofIndex = 0; dofIndex < nDofsPerElement; dofIndex++)
    {
      for (int dofComponent = 0; dofComponent < D; dofComponent++)
      {
        // compute index of degree of freedom and component (integrade values vector index)
        const int i = D*dofIndex + dofComponent;

        // get integrated value
        double integratedValue = integratedValues[i];

        LOG(DEBUG) << "  dof " << dofIndex << ", component " << dofComponent << " integrated value: " << integratedValue;

        // Compute indices in result vector. For each dof there are D values for the displacement components.
        // Therefore the nDofsPerElement number is not the number of unknows.
        dof_no_t resultVectorIndex = dofNo[dofIndex]*D + dofComponent;

        ierr = VecSetValue(resultVec, resultVectorIndex, integratedValue, ADD_VALUES); CHKERRV(ierr);
        //LOG(DEBUG) << integratedValue << ", " << resultVectorIndex;

      }  // dofComponent
    }  // dofIndex
  }  // elementNo

  // -------------- surface traction in current configuration --------------
  for (typename std::vector<typename SolidMechanicsBoundaryConditions<BasisOnMeshType,Term>::TractionBoundaryCondition>::const_iterator iter = this->tractionCurrentConfiguration_.begin();
       iter != this->tractionCurrentConfiguration_.end(); iter++)
  {
    /*
    element_no_t elementGlobalNo;

    Mesh::face_t face;
    std::vector<std::pair<dof_no_t, Vec3>> dofVectors;  //<element-local dof no, value>
    */
    element_no_t elementNo = iter->elementGlobalNo;

    // check if element no is valid
    if (elementNo < 0 || elementNo > mesh->nElements())
    {
      LOG(ERROR) << "Element " << elementNo << " for which surface traction (cur.conf.) is specified is invalid (number of elements: " << mesh->nElements() << ")";
      continue;
    }

    // get geometry field of reference configuration, this is always 3D also for 2D problems
    std::array<Vec3,nDofsPerElement> geometryReferenceValues;
    this->data_.geometryReference().getElementValues(elementNo, geometryReferenceValues);

    // get displacement field
    std::array<VecD<D>,nDofsPerElement> displacementValues;
    this->data_.displacements().getElementValues(elementNo, displacementValues);

    // loop over integration points (e.g. gauss points)
    for (unsigned int samplingPointIndex = 0; samplingPointIndex < samplingPointsSurface.size(); samplingPointIndex++)
    {
      // get parameter values of current sampling point
      std::array<double,D-1> xiSurface = samplingPointsSurface[samplingPointIndex];
      VecD<D> xi = Mesh::getXiOnFace(iter->face, xiSurface);

      // compute the normal in reference configuration of the specified face at point xi
      VecD<D> normalReferenceConfiguration = MathUtility::transformToD<D,3>(mesh->getNormal(iter->face, geometryReferenceValues, xi));

      // compute the 3xD jacobian of the parameter space to world space mapping
      Tensor2<D> jacobianMaterial = MathUtility::transformToDxD<D,D>(BasisOnMesh::computeJacobian(geometryReferenceValues, xi));
      double jacobianDeterminant;
      Tensor2<D> inverseJacobianMaterial = MathUtility::computeInverse<D>(jacobianMaterial, jacobianDeterminant);
      // jacobianMaterial[columnIdx][rowIdx] = dX_rowIdx/dxi_columnIdx
      // inverseJacobianMaterial[columnIdx][rowIdx] = dxi_rowIdx/dX_columnIdx because of inverse function theorem

      // F
      Tensor2<D> deformationGradient = this->computeDeformationGradient(displacementValues, inverseJacobianMaterial, xi);

      // J*F^{-T}
      Tensor2<D> cofactorMatrix = MathUtility::computeCofactorMatrix<D>(deformationGradient);

      VecD<D> normalCurrentConfiguration = cofactorMatrix * normalReferenceConfiguration;
      double surfaceStretch = MathUtility::norm<D>(normalCurrentConfiguration);  // ds/dS, because |dS| = |normalReferenceConfiguration| = 1

      LOG(DEBUG) << "surfaceStretch: " << surfaceStretch;

      // T = t ds/dS (ds, dS are scalar)
      // Nansons formula: ds = J F^-T dS (ds, dS are vectors), ds/dS = |ds|/|dS|  = |J F^-T dS|/|dS|, construct |dS| = 1

      // set all entries to 0
      evaluationsArraySurface[samplingPointIndex] = {0.0};

      // compute the traction at xi by summing contributions from elemental dofs
      VecD<D> tractionCurrentConfiguration({0.0});
      for (typename std::vector<std::pair<dof_no_t, VecD<D>>>::const_iterator dofVectorsIter = iter->dofVectors.begin();
           dofVectorsIter != iter->dofVectors.end();
           dofVectorsIter++)
      {
        int dofIndex = dofVectorsIter->first;
        VecD<D> forceVector = dofVectorsIter->second;

        tractionCurrentConfiguration += forceVector * mesh->phi(dofIndex, xi);
      }

      // loop over dofs of element
      for (typename std::vector<std::pair<dof_no_t, VecD<D>>>::const_iterator dofVectorsIter = iter->dofVectors.begin();
           dofVectorsIter != iter->dofVectors.end();
           dofVectorsIter++)
      {
        int dofIndex = dofVectorsIter->first;

        VecD<D> dofIntegrand = tractionCurrentConfiguration * surfaceStretch * MathUtility::sqr(mesh->phi(dofIndex, xi));

        // store integrand in evaluations array
        for (int i = 0; i < D; i++)
        {
          evaluationsArraySurface[samplingPointIndex][dofIndex*D+i] = dofIntegrand[i];
        }
      }  // dofIndex

    }  // function evaluations

    // integrate all values for result vector entries at once
    EvaluationsType integratedValues = QuadratureSurface::computeIntegral(evaluationsArraySurface);

    // get indices of element-local dofs
    std::array<dof_no_t,nDofsPerElement> dofNo = mesh->getElementDofNos(elementNo);

    // add entries in result vector
    // loop over indices of unknows (dofIndex,dofComponent)
    for (int dofIndex = 0; dofIndex < nDofsPerElement; dofIndex++)
    {
      for (int dofComponent = 0; dofComponent < D; dofComponent++)
      {
        // compute index of degree of freedom and component (integrade values vector index)
        const int i = D*dofIndex + dofComponent;

        // get integrated value
        double integratedValue = integratedValues[i];

        // Compute indices in result vector. For each dof there are D values for the displacement components.
        // Therefore the nDofsPerElement number is not the number of unknows.
        dof_no_t resultVectorIndex = dofNo[dofIndex]*D + dofComponent;

        ierr = VecSetValue(resultVec, resultVectorIndex, integratedValue, ADD_VALUES); CHKERRV(ierr);

      }  // dofComponent
    }  // dofIndex
  }  // elementNo

  ierr = VecAssemblyBegin(resultVec); CHKERRV(ierr);
  ierr = VecAssemblyEnd(resultVec); CHKERRV(ierr);
  // return memory acces of result vector back to PETSc (not used)
  //VecRestoreArray(resultVec, &result);

  if (VLOG_IS_ON(1)) VLOG(1) << "  ->wExt: " << PetscUtility::getStringVector(resultVec);
}

template<typename BasisOnMeshType,typename BasisOnMeshTypeForUtility, typename MixedQuadratureType, typename Term>
void SolidMechanicsCommon<BasisOnMeshType,BasisOnMeshTypeForUtility,MixedQuadratureType,Term>::
computeInternalVirtualWork(Vec &resultVec)
{
  LOG(TRACE) << "computeInternalVirtualWork";

  typedef typename BasisOnMeshType::HighOrderBasisOnMesh BasisOnMesh;  // for mixed formulation get the high order BasisOnMesh
  typedef typename MixedQuadratureType::HighOrderQuadrature QuadratureType;  // for mixed formulation get the high order quadrature

  // get pointer to mesh object
  std::shared_ptr<BasisOnMesh> mesh = this->data_.mesh();

  const int D = BasisOnMesh::dim();  // = 2 or 3
  //const int nDofs = mesh->nDofs();
  const int nDofsPerElement = BasisOnMesh::nDofsPerElement();
  const int nElements = mesh->nElements();
  const int nUnknowsPerElement = nDofsPerElement*D;    // D directions for displacements per dof

  // define shortcuts for quadrature
  typedef Quadrature::TensorProduct<D,QuadratureType> QuadratureDD;

  // define type to hold evaluations of integrand
  typedef std::array<double, nUnknowsPerElement> EvaluationsType;
  typedef std::array<
            EvaluationsType,
            QuadratureDD::numberEvaluations()
          > EvaluationsArrayType;     // evaluations[nGP^D](nUnknows)

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

  // set values to zero
  ierr = VecZeroEntries(resultVec); CHKERRV(ierr);

  // get memory from PETSc where to store result (not used)
  //const PetscScalar *result;
  //VecGetArray(resultVec, &result);

  // loop over elements
  for (int elementNo = 0; elementNo < nElements; elementNo++)
  {

#ifdef QUADRATURE_TEST
    Control::PerformanceMeasurement::start("internalVirtualWork");
#endif

    // get geometry field of reference configuration, note the dimension of the vecs is always 3 also for 2D problems
    std::array<Vec3,nDofsPerElement> geometryReferenceValues;
    this->data_.geometryReference().getElementValues(elementNo, geometryReferenceValues);

    // get displacement field values for element
    std::array<VecD<D>,nDofsPerElement> displacementValues;
    this->data_.displacements().getElementValues(elementNo, displacementValues);

    // For mixed formulation get the pressure values of this element. This is done only in the derived class for mixed formulation and not in the derived class for penalty formulation.
    this->preparePressureInterpolation(elementNo);

    // loop over integration points (e.g. gauss points) for displacement field
    for (unsigned int samplingPointIndex = 0; samplingPointIndex < samplingPoints.size(); samplingPointIndex++)
    {
      // get parameter values of current sampling point
      VecD<D> xi = samplingPoints[samplingPointIndex];

      // compute the 3xD jacobian of the parameter space to world space mapping
      Tensor2<D> jacobianMaterial = MathUtility::transformToDxD<D,D>(BasisOnMesh::computeJacobian(geometryReferenceValues, xi));
      double jacobianDeterminant;
      Tensor2<D> inverseJacobianMaterial = MathUtility::computeInverse<D>(jacobianMaterial, jacobianDeterminant);
      // jacobianMaterial[columnIdx][rowIdx] = dX_rowIdx/dxi_columnIdx
      // inverseJacobianMaterial[columnIdx][rowIdx] = dxi_rowIdx/dX_columnIdx because of inverse function theorem

      checkInverseIsCorrect<D>(jacobianMaterial, inverseJacobianMaterial, "jacobianMaterial");

      // F
      Tensor2<D> deformationGradient = this->computeDeformationGradient(displacementValues, inverseJacobianMaterial, xi);
      double deformationGradientDeterminant = MathUtility::computeDeterminant<D>(deformationGradient);  // J

      Tensor2<D> rightCauchyGreen = this->computeRightCauchyGreenTensor(deformationGradient);  // C = F^T*F

      checkSymmetry<D>(rightCauchyGreen, "C");

      double rightCauchyGreenDeterminant;   // J^2
      Tensor2<D> inverseRightCauchyGreen = MathUtility::computeSymmetricInverse<D>(rightCauchyGreen, rightCauchyGreenDeterminant);  // C^-1

      checkInverseIsCorrect<D>(rightCauchyGreen, inverseRightCauchyGreen, "rightCauchyGreen");

      // invariants
      std::array<double,3> invariants = this->computeInvariants(rightCauchyGreen, rightCauchyGreenDeterminant);  // I_1, I_2, I_3
      std::array<double,2> reducedInvariants = this->computeReducedInvariants(invariants, deformationGradientDeterminant); // Ibar_1, Ibar_2

      // pressure is the artificialPressure for penalty formulation or the separately interpolated pressure for mixed formulation
      double pressureTilde = 0.0;
      const double pressure = this->getPressure(deformationGradientDeterminant, xi, pressureTilde);

      // Pk2 stress tensor S = S_vol + S_iso (p.234)
      Tensor2<D> fictitiousPK2Stress;
      Tensor2<D> pk2StressIsochoric;

      //! compute 2nd Piola-Kirchhoff stress tensor S = 2*dPsi/dC and the fictitious PK2 Stress Sbar
      Tensor2<D> PK2Stress = this->computePK2Stress(pressure, rightCauchyGreen, inverseRightCauchyGreen, reducedInvariants, deformationGradientDeterminant, fictitiousPK2Stress, pk2StressIsochoric);

      // check if PK2 is the same as from explicit formula for Mooney-Rivlin (p.249)
      if (VLOG_IS_ON(1))
      {
        this->checkFictitiousPK2Stress(fictitiousPK2Stress, rightCauchyGreen, deformationGradientDeterminant, reducedInvariants);
      }

      std::array<VecD<D>,nDofsPerElement> gradPhi = mesh->getGradPhi(xi);
      // (column-major storage) gradPhi[M][a] = dphi_M / dxi_a
      // gradPhi[column][row] = gradPhi[dofIndex][i] = dphi_dofIndex/dxi_i, columnIdx = dofIndex, rowIdx = which direction

      VLOG(2) << "";
      VLOG(2) << "element " << elementNo << " xi: " << xi;
      VLOG(2) << "  geometryReferenceValues: " << geometryReferenceValues;
      VLOG(2) << "  displacementValues: " << displacementValues;
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
      VLOG(2) << "  fictitiousPK2Stress: Sbar=" << fictitiousPK2Stress;
      VLOG(2) << "  pk2StressIsochoric: S_iso=" << pk2StressIsochoric;
      VLOG(2) << "  PK2Stress: S=" << PK2Stress;
      VLOG(2) << "  gradPhi: " << gradPhi;

      VLOG(1) << "  xi: " << xi << ", J: " << deformationGradientDeterminant << ", p: " << pressure << ", S11: " << PK2Stress[0][0];

      if (samplingPointIndex == 0 && D == 3)
        LOG(DEBUG) << " F11: " << deformationGradient[0][0] << ", F22,F33: " << deformationGradient[1][1] <<"," << deformationGradient[2][2]
          << ", F12,F13,F23: " << deformationGradient[1][0] << "," << deformationGradient[2][0] << "," << deformationGradient[2][1]
          << ", J: " << deformationGradientDeterminant << ", p: " << pressure;

      if (VLOG_IS_ON(2))
      {
        Tensor2<D> greenLangrangeStrain = this->computeGreenLagrangeStrain(rightCauchyGreen);
        VLOG(2) << "  strain E=" << greenLangrangeStrain;
      }

      if (deformationGradientDeterminant < 1e-12)   // if any entry of the deformation gradient is negative
      {
        LOG(FATAL) << "Deformation gradient " << deformationGradient << " has zero or negative determinant " << deformationGradientDeterminant;
      }

      // loop basis functions and evaluate integrand at xi
      for (int aDof = 0; aDof < nDofsPerElement; aDof++)           // index over dofs, each dof has D components, L in derivation
      {
        for (int aComponent = 0; aComponent < D; aComponent++)     // lower-case a in derivation, index over displacement components
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
                const double dxik_dXA = inverseJacobianMaterial[aInternal][k];  // inverseJacobianMaterial[A][k] = J_kA = dxi_k/dX_A

                dphiL_dXA += dphiL_dxik * dxik_dXA;

                // compute dphiL/dXB from dphiL/dxik and dxik/dXB
                const double dxik_dXB = inverseJacobianMaterial[bInternal][k];  // inverseJacobianMaterial[B][k] = J_kB = dxi_k/dX_B

                dphiL_dXB += dphiL_dxik * dxik_dXB;
              }

              integrand += 1./2. * PK2Stress[bInternal][aInternal] * (faB*dphiL_dXA + faA*dphiL_dXB);

            }  // B, bInternal
          }  // A, aInternal

          VLOG(2) << "   (L,a)=(" << aDof << "," << aComponent << "), integrand: " << integrand;

          // store integrand in evaluations array
          evaluationsArray[samplingPointIndex][i] = integrand;

        }  // a
      }  // L
    }  // function evaluations

    // integrate all values for result vector entries at once
    EvaluationsType integratedValues = QuadratureDD::computeIntegral(evaluationsArray);

#ifdef QUADRATURE_TEST
    Control::PerformanceMeasurement::stop("internalVirtualWork");
#endif


#ifdef QUADRATURE_TEST

    // loop over integration points (e.g. gauss points) for displacement field
    for (unsigned int samplingPointIndex = 0; samplingPointIndex < samplingPointsExact.size(); samplingPointIndex++)
    {
      // get parameter values of current sampling point
      VecD<D> xi = samplingPointsExact[samplingPointIndex];

      // compute the 3xD jacobian of the parameter space to world space mapping
      Tensor2<D> jacobianMaterial = MathUtility::transformToDxD<D,D>(BasisOnMesh::computeJacobian(geometryReferenceValues, xi));
      double jacobianDeterminant;
      Tensor2<D> inverseJacobianMaterial = MathUtility::computeInverse<D>(jacobianMaterial, jacobianDeterminant);
      // jacobianMaterial[columnIdx][rowIdx] = dX_rowIdx/dxi_columnIdx
      // inverseJacobianMaterial[columnIdx][rowIdx] = dxi_rowIdx/dX_columnIdx because of inverse function theorem

      checkInverseIsCorrect<D>(jacobianMaterial, inverseJacobianMaterial, "jacobianMaterial");

      // F
      Tensor2<D> deformationGradient = this->computeDeformationGradient(displacementValues, inverseJacobianMaterial, xi);
      double deformationGradientDeterminant = MathUtility::computeDeterminant<D>(deformationGradient);  // J

      Tensor2<D> rightCauchyGreen = this->computeRightCauchyGreenTensor(deformationGradient);  // C = F^T*F

      checkSymmetry<D>(rightCauchyGreen, "C");

      double rightCauchyGreenDeterminant;   // J^2
      Tensor2<D> inverseRightCauchyGreen = MathUtility::computeSymmetricInverse<D>(rightCauchyGreen, rightCauchyGreenDeterminant);  // C^-1

      checkInverseIsCorrect<D>(rightCauchyGreen, inverseRightCauchyGreen, "rightCauchyGreen");

      // invariants
      std::array<double,3> invariants = this->computeInvariants(rightCauchyGreen, rightCauchyGreenDeterminant);  // I_1, I_2, I_3
      std::array<double,2> reducedInvariants = this->computeReducedInvariants(invariants, deformationGradientDeterminant); // Ibar_1, Ibar_2

      // pressure is the artificialPressure for penalty formulation or the separately interpolated pressure for mixed formulation
      double pressureTilde = 0.0;
      const double pressure = this->getPressure(deformationGradientDeterminant, xi, pressureTilde);

      // Pk2 stress tensor S = S_vol + S_iso (p.234)
      Tensor2<D> fictitiousPK2Stress;
      Tensor2<D> pk2StressIsochoric;

      //! compute 2nd Piola-Kirchhoff stress tensor S = 2*dPsi/dC and the fictitious PK2 Stress Sbar
      Tensor2<D> PK2Stress = this->computePK2Stress(pressure, rightCauchyGreen, inverseRightCauchyGreen, reducedInvariants, deformationGradientDeterminant, fictitiousPK2Stress, pk2StressIsochoric);

      // check if PK2 is the same as from explicit formula for Mooney-Rivlin (p.249)
      if (VLOG_IS_ON(1))
      {
        this->checkFictitiousPK2Stress(fictitiousPK2Stress, rightCauchyGreen, deformationGradientDeterminant, reducedInvariants);
      }

      std::array<VecD<D>,nDofsPerElement> gradPhi = mesh->getGradPhi(xi);
      // (column-major storage) gradPhi[M][a] = dphi_M / dxi_a
      // gradPhi[column][row] = gradPhi[dofIndex][i] = dphi_dofIndex/dxi_i, columnIdx = dofIndex, rowIdx = which direction

      // loop basis functions and evaluate integrand at xi
      for (int aDof = 0; aDof < nDofsPerElement; aDof++)           // index over dofs, each dof has D components, L in derivation
      {
        for (int aComponent = 0; aComponent < D; aComponent++)     // lower-case a in derivation, index over displacement components
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
                const double dxik_dXA = inverseJacobianMaterial[aInternal][k];  // inverseJacobianMaterial[A][k] = J_kA = dxi_k/dX_A

                dphiL_dXA += dphiL_dxik * dxik_dXA;

                // compute dphiL/dXB from dphiL/dxik and dxik/dXB
                const double dxik_dXB = inverseJacobianMaterial[bInternal][k];  // inverseJacobianMaterial[B][k] = J_kB = dxi_k/dX_B

                dphiL_dXB += dphiL_dxik * dxik_dXB;
              }

              integrand += 1./2. * PK2Stress[bInternal][aInternal] * (faB*dphiL_dXA + faA*dphiL_dXB);

            }  // B, bInternal
          }  // A, aInternal

          // store integrand in evaluations array
          evaluationsExactArray[samplingPointIndex][i] = integrand;

        }  // a
      }  // L
    }  // function evaluations

    // integrate all values for result vector entries at once
    EvaluationsType integratedValuesExact = QuadratureExactDD::computeIntegral(evaluationsExactArray);

    Control::PerformanceMeasurement::measureError("internalVirtualWork", (integratedValues - integratedValuesExact)/integratedValuesExact);

#endif

    // get indices of element-local dofs
    std::array<dof_no_t,nDofsPerElement> dofNo = mesh->getElementDofNos(elementNo);

    VLOG(2) << "  element " << elementNo << " has dofs " << dofNo;

    // add entries in result vector
    // loop over indices of unknows (aDof,aComponent)
    for (int aDof = 0; aDof < nDofsPerElement; aDof++)
    {
      for (int aComponent = 0; aComponent < D; aComponent++)
      {
        // compute index of degree of freedom and component (matrix row index)
        const int i = D*aDof + aComponent;

        // integrate value and set entry in stiffness matrix
        double integratedValue = integratedValues[i];

        // Compute indices in stiffness matrix. For each dof there are D values for the D displacement components
        // Therefore the nDofsPerElement number is not the number of unknows.
        dof_no_t resultVectorIndex = dofNo[aDof]*D + aComponent;

        VLOG(2) << "  result vector (L,a)=("<<aDof<<","<<aComponent<<"), " <<i<<", dof " << dofNo[aDof];
        VLOG(2) << "      vector index (unknown no): "<<resultVectorIndex<<", integrated value: "<<integratedValue;

        ierr = VecSetValue(resultVec, resultVectorIndex, integratedValue, ADD_VALUES); CHKERRV(ierr);

      }  // aComponent
    }  // aDof
  }  // elementNo

  ierr = VecAssemblyBegin(resultVec); CHKERRV(ierr);
  ierr = VecAssemblyEnd(resultVec); CHKERRV(ierr);
  // return memory acces of result vector back to PETSc (not used)
  //VecRestoreArray(resultVec, &result);

  //if (VLOG_IS_ON(1)) VLOG(1) << "    disp: " << PetscUtility::getStringVector(this->data_.displacements().values());
  if (VLOG_IS_ON(1)) VLOG(1) << "  ->wInt: " << PetscUtility::getStringVector(resultVec);
}


template<typename BasisOnMeshType,typename BasisOnMeshTypeForUtility, typename QuadratureType, typename Term>
void SolidMechanicsCommon<BasisOnMeshType,BasisOnMeshTypeForUtility,QuadratureType,Term>::
getExternalVirtualWork(Vec &resultVec)
{
  if (this->data_.externalVirtualWorkIsConstant())
  {
    // the external virtual work is independent of the current displacements and was already computed from initializeBoundaryConditions
    VecCopy(this->data_.externalVirtualWork().values(), resultVec);
  }
  else
  {
    computeExternalVirtualWork(resultVec);
    this->data_.externalVirtualWork().values() = resultVec;
  }

}

template<typename BasisOnMeshType,typename BasisOnMeshTypeForUtility, typename QuadratureType, typename Term>
void SolidMechanicsCommon<BasisOnMeshType,BasisOnMeshTypeForUtility,QuadratureType,Term>::
computeInternalMinusExternalVirtualWork(Vec &resultVec)
{
  PetscErrorCode ierr;

  if (this->data_.computeWithReducedVectors())
  {
    Vec &wInt = this->data_.internalVirtualWork().values();
    Vec &wExt = this->data_.externalVirtualWork().values();
    Vec &wIntReduced = this->data_.internalVirtualWorkReduced();

    const int D = BasisOnMeshType::dim();
    const int nUnknowns = this->data_.mesh()->nDofs() * D;

    this->computeInternalVirtualWork(wInt);
    LOG(DEBUG) << "--                           dW_int: " << PetscUtility::getStringVector(wInt);

    this->getExternalVirtualWork(wExt);

    // remove entries that correspond to Dirichlet BC
    this->reduceVector(wInt, wIntReduced, nUnknowns);
    this->reduceVector(wExt, resultVec, nUnknowns);

    LOG(DEBUG) << "--                           dW_int: " << PetscUtility::getStringVector(wInt);
    LOG(DEBUG) << "--                           dW_ext: " << PetscUtility::getStringVector(wExt);

    // compute wInt - wExt
    ierr = VecAYPX(resultVec, -1, wIntReduced); CHKERRV(ierr);
  }
  else
  {
    Vec &wInt = this->data_.internalVirtualWork().values();

    computeInternalVirtualWork(wInt);
    getExternalVirtualWork(resultVec);

    LOG(DEBUG) << "--                           dW_int: " << PetscUtility::getStringVector(wInt);
    LOG(DEBUG) << "--                           dW_ext: " << PetscUtility::getStringVector(resultVec);

    // compute wInt - wExt
    ierr = VecAYPX(resultVec, -1, wInt); CHKERRV(ierr);
  }
}

template<typename BasisOnMeshType,typename BasisOnMeshTypeForUtility, typename MixedQuadratureType, typename Term>
void SolidMechanicsCommon<BasisOnMeshType,BasisOnMeshTypeForUtility,MixedQuadratureType,Term>::
initialize()
{
  initializeMaterialParameters();
  this->data_.initialize();

  // initialize external virtual energy
  this->initializeBoundaryConditions(this->data_.externalVirtualWorkIsConstant(), this->data_.nUnknowns(), this->specificSettings_, this->data_);

  // nUnknowns gives the number of unknowns including Dirichlet BC values, for mixed formulation number of displacement and pressure unknowns, for penalty formulation only displacements unknowns
  const int vectorSize = this->nUnknowns() - this->dirichletValues_.size();

  // initialize Vecs that are used by nonlinear solver
  this->data_.initializeSolverVariables(vectorSize);

  // if the external virtual work is constant, i.e. independent of the displacement, it does not have to be computed for every new iteration of the nonlinear solver. Compute it once, now.
  if (this->data_.externalVirtualWorkIsConstant())
  {
    this->computeExternalVirtualWork(this->data_.externalVirtualWork().values());

    VLOG(1) << "-- initially computed dW_ext: " << PetscUtility::getStringVector(this->data_.externalVirtualWork().values());
  }

  this->outputIntermediateSteps_ = PythonUtility::getOptionBool(this->specificSettings_, "outputIntermediateSteps", false);

  this->printBoundaryConditions();

  //this->setStiffnessMatrix();
  //this->setRightHandSide();
  this->data_.finalAssembly();
}

template<typename BasisOnMeshType,typename BasisOnMeshTypeForUtility, typename QuadratureType, typename Term>
void SolidMechanicsCommon<BasisOnMeshType,BasisOnMeshTypeForUtility,QuadratureType,Term>::
initializeMaterialParameters()
{
  std::array<double,Term::nMaterialParameters> parameters
    = PythonUtility::template getOptionArray<double, Term::nMaterialParameters>(this->specificSettings_, "materialParameters", 0.0);

  LOG(DEBUG) << "Material has " << Term::nMaterialParameters << " parameters, parsed from \"materialParameters\": " << parameters;

  std::vector<double> parametersVector(Term::nMaterialParameters);
  std::copy(parameters.begin(), parameters.end(), parametersVector.begin());

  // set all PARAM(i) values to the values given by materialParameters
  SEMT::set_parameters<Term::nMaterialParameters>::to(parametersVector);
}

template<typename BasisOnMeshType,typename BasisOnMeshTypeForUtility, typename QuadratureType, typename Term>
bool SolidMechanicsCommon<BasisOnMeshType,BasisOnMeshTypeForUtility,QuadratureType,Term>::
computeWithReducedVectors()
{
  return this->data_.computeWithReducedVectors();
}

template<typename BasisOnMeshType,typename BasisOnMeshTypeForUtility, typename QuadratureType, typename Term>
void SolidMechanicsCommon<BasisOnMeshType,BasisOnMeshTypeForUtility,QuadratureType,Term>::
computeAnalyticStiffnessMatrix(Mat &solverStiffnessMatrix)
{
  if (this->data_.computeWithReducedVectors())
  {
    // for mixed formulation sum of number of u and p unknowns
    const int nUnknownsFull = this->nUnknowns();

    // Compute new stiffness matrix from current displacements in this->data_.tangentStiffnessMatrix(). This uses the full vectors.
    this->setStiffnessMatrix();

    // create new reduced stiffness matrix by copying the non-zero entries
    this->reduceMatrix(this->data_.tangentStiffnessMatrix(), solverStiffnessMatrix, nUnknownsFull);
  }
  else
  {
    // Compute new stiffness matrix from current displacements. This uses the full vectors.
    this->setStiffnessMatrix(solverStiffnessMatrix);

    // set the appropriate entries to match Dirichlet BC
    this->applyDirichletBoundaryConditionsInStiffnessMatrix(solverStiffnessMatrix, this->data_);
  }
}

template<typename BasisOnMeshType,typename BasisOnMeshTypeForUtility, typename QuadratureType, typename Term>
void SolidMechanicsCommon<BasisOnMeshType,BasisOnMeshTypeForUtility,QuadratureType,Term>::
writeOutput()
{
  if (this->outputIntermediateSteps_)
  {
    this->updateGeometryActual();

    static int fileCounter = 0;
    this->outputWriterManager_.writeOutput(this->data_, fileCounter);
    fileCounter++;
  }

#ifdef QUADRATURE_TEST
  Control::PerformanceMeasurement::log("timing.txt");
#endif
}

};    // namespace