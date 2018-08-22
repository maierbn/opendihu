#include "spatial_discretization/finite_element_method/03_assemble_rhs.h"

#include <Python.h>
#include <memory>
#include <vector>
#include <petscsys.h>

#include "quadrature/tensor_product.h"
#include "basis_on_mesh/basis_on_mesh.h"
#include "spatial_discretization/finite_element_method/03_integrand_rhs.h"
#include "field_variable/field_variable.h"

namespace SpatialDiscretization
{

// 1D,2D,3D rhs vector of Deformable mesh
template<typename BasisOnMeshType, typename QuadratureType, typename Term, typename Dummy>
void AssembleRightHandSide<BasisOnMeshType, QuadratureType, Term, Dummy>::
transferRhsToWeakForm()
{
  const int D = BasisOnMeshType::dim();
  LOG(TRACE)<<"transferRhsToWeakForm " << D << "D";

  // define shortcuts for integrator and basis
  typedef Quadrature::TensorProduct<D,QuadratureType> QuadratureDD;
  const int nDofsPerElement = BasisOnMeshType::nDofsPerElement();
  typedef MathUtility::Matrix<nDofsPerElement,nDofsPerElement> EvaluationsType;
  typedef std::array<
            EvaluationsType,
            QuadratureDD::numberEvaluations()
          > EvaluationsArrayType;    // evaluations[nGP^D][nDofs][nDofs]

  // initialize variables
  FieldVariable::FieldVariable<BasisOnMeshType,1> &rightHandSide = this->data_.rightHandSide();

  std::shared_ptr<BasisOnMeshType> mesh = std::static_pointer_cast<BasisOnMeshType>(this->data_.mesh());

  // merge local changes on the vector
  rightHandSide.startVectorManipulation();
  
  // get all entries
  std::vector<double> rhsValues;
  rightHandSide.getLocalValues(rhsValues);

  // initialize values to zero
  rightHandSide.zeroEntries();

  // setup arrays used for integration
  std::array<std::array<double,D>, QuadratureDD::numberEvaluations()> samplingPoints = QuadratureDD::samplingPoints();
  EvaluationsArrayType evaluationsArray{};

  LOG(DEBUG) << "1D integration with " << QuadratureType::numberEvaluations() << " evaluations";
  LOG(DEBUG) << D << "D integration with " << QuadratureDD::numberEvaluations() << " evaluations";

  // set entries in rhs vector
  // loop over elements
  for (element_no_t elementNo = 0; elementNo < mesh->nElementsLocal(); elementNo++)
  {
    // get indices of element-local dofs
    std::array<dof_no_t,nDofsPerElement> dofNosLocal = mesh->getElementDofNosLocal(elementNo);

    VLOG(2) << "element " << elementNo;

    // get geometry field (which are the node positions for Lagrange basis and node positions and derivatives for Hermite)
    std::array<Vec3,BasisOnMeshType::nDofsPerElement()> geometry;
    mesh->getElementGeometry(elementNo, geometry);

    // compute integral
    for (unsigned int samplingPointIndex = 0; samplingPointIndex < samplingPoints.size(); samplingPointIndex++)
    {
      // evaluate function to integrate at samplingPoints[i*2], write value to evaluations[i]
      std::array<double,D> xi = samplingPoints[samplingPointIndex];

      // compute the 3xD jacobian of the parameter space to world space mapping
      auto jacobian = BasisOnMeshType::computeJacobian(geometry, xi);

      // get evaluations of integrand which is defined in another class
      evaluationsArray[samplingPointIndex] = IntegrandRightHandSide<D,EvaluationsType,BasisOnMeshType,Term>::evaluateIntegrand(jacobian,xi);

    }  // function evaluations

    // integrate all values for the (i,j) dof pairs at once
    EvaluationsType integratedValues = QuadratureDD::computeIntegral(evaluationsArray);

    // perform integration and add to entry in rhs vector
    for (int i=0; i<nDofsPerElement; i++)
    {
      for (int j=0; j<nDofsPerElement; j++)
      {
        // integrate value and set entry in stiffness matrix
        double integratedValue = integratedValues(i,j);

        double value = integratedValue * rhsValues[dofNosLocal[j]];
        VLOG(2) << "  dof pair (" << i<<","<<j<<"), integrated value: "<<integratedValue<<", rhsValue["<<dofNosLocal[j]<<"]: " << rhsValues[dofNosLocal[j]] <<" = " << value;

        rightHandSide.setValue(dofNosLocal[i], value, ADD_VALUES);
      }  // j
    }  // i
  }  // elementNo

  // merge local changes on the vector, parallel assembly
  rightHandSide.finishVectorManipulation();
 
}

// 1D,2D,3D rhs discretization matrix, i.e. matrix that transforms rhs values to discretized form, of Deformable mesh
template<typename BasisOnMeshType, typename QuadratureType, typename Term, typename Dummy>
void AssembleRightHandSide<BasisOnMeshType, QuadratureType, Term, Dummy>::
setMassMatrix()
{
  // check if matrix discretization matrix exists
  if (!this->data_.massMatrixInitialized())
  {
    this->data_.initializeMassMatrix();

    const int D = BasisOnMeshType::dim();
    LOG(TRACE)<<"createMassMatrix " << D << "D";

    // massMatrix * f_strong = rhs_weak
    // row of massMatrix: contributions to a single entry in rhs_weak

    // define shortcuts for integrator and basis
    typedef Quadrature::TensorProduct<D,QuadratureType> QuadratureDD;
    const int nDofsPerElement = BasisOnMeshType::nDofsPerElement();
    typedef MathUtility::Matrix<nDofsPerElement,nDofsPerElement> EvaluationsType;
    typedef std::array<
              EvaluationsType,
              QuadratureDD::numberEvaluations()
            > EvaluationsArrayType;     // evaluations[nGP^D][nDofs][nDofs]

    // initialize variables
    std::shared_ptr<PartitionedPetscMat<BasisOnMeshType>> massMatrix = this->data_.massMatrix();

    std::shared_ptr<BasisOnMeshType> mesh = std::static_pointer_cast<BasisOnMeshType>(this->data_.mesh());

    // initialize values to zero
    // loop over elements
    for (element_no_t elementNo = 0; elementNo < mesh->nElementsLocal(); elementNo++)
    {
      std::array<dof_no_t,nDofsPerElement> dofNosLocal = mesh->getElementDofNosLocal(elementNo);

      for (int i=0; i<nDofsPerElement; i++)
      {
        for (int j=0; j<nDofsPerElement; j++)
        {
          massMatrix->setValue(dofNosLocal[i], dofNosLocal[j], 0, INSERT_VALUES);
        }
      }
    }

    // setup arrays used for integration
    std::array<std::array<double,D>, QuadratureDD::numberEvaluations()> samplingPoints = QuadratureDD::samplingPoints();
    EvaluationsArrayType evaluationsArray{};

    LOG(DEBUG) << "1D integration with " << QuadratureType::numberEvaluations() << " evaluations";
    LOG(DEBUG) << D << "D integration with " << QuadratureDD::numberEvaluations() << " evaluations";
    
    // set entries in massMatrix
    // loop over elements
    for (element_no_t elementNo = 0; elementNo < mesh->nElementsLocal(); elementNo++)
    {
      // get indices of element-local dofs
      std::array<dof_no_t,nDofsPerElement> dofNosLocal = mesh->getElementDofNosLocal(elementNo);

      // get geometry field (which are the node positions for Lagrange basis and node positions and derivatives for Hermite)
      std::array<Vec3,BasisOnMeshType::nDofsPerElement()> geometry;
      mesh->getElementGeometry(elementNo, geometry);

      // compute integral
      for (unsigned int samplingPointIndex = 0; samplingPointIndex < samplingPoints.size(); samplingPointIndex++)
      {
        // evaluate function to integrate at samplingPoints[i*2], write value to evaluations[i]
        std::array<double,D> xi = samplingPoints[samplingPointIndex];

        // compute the 3xD jacobian of the parameter space to world space mapping
        auto jacobian = BasisOnMeshType::computeJacobian(geometry, xi);

        // get evaluations of integrand which is defined in another class
        evaluationsArray[samplingPointIndex] = IntegrandRightHandSide<D,EvaluationsType,BasisOnMeshType,Term>::evaluateIntegrand(jacobian,xi);

      }  // function evaluations

      // integrate all values for the (i,j) dof pairs at once
      EvaluationsType integratedValues = QuadratureDD::computeIntegral(evaluationsArray);

      // perform integration and add to entry in rhs vector
      for (int i=0; i<nDofsPerElement; i++)
      {
        for (int j=0; j<nDofsPerElement; j++)
        {
          // integrate value and set entry in discretization matrix
          double integratedValue = integratedValues(i,j);

          massMatrix->setValue(dofNosLocal[i], dofNosLocal[j], integratedValue, ADD_VALUES);
        }  // j
      }  // i
    }  // elementNo

    // merge local changes in parallel and assemble the matrix (MatAssemblyBegin, MatAssemblyEnd)
    massMatrix->assembly(MAT_FINAL_ASSEMBLY);
  }
}

/*
template<typename MixedBasisOnMeshType, typename MixedQuadratureType>
AssembleRightHandSide<MixedBasisOnMeshType, MixedQuadratureType, Equation::Static::SolidMechanics>::
AssembleRightHandSide(DihuContext &context) :
  FiniteElementMethodStiffnessMatrix<MixedBasisOnMeshType, MixedQuadratureType, Equation::Static::SolidMechanics>::FiniteElementMethodStiffnessMatrix(context)
{

}*/

};