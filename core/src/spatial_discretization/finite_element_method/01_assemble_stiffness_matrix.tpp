#include "spatial_discretization/finite_element_method/01_assemble_stiffness_matrix.h"

#include <Python.h>
#include <memory>
#include <vector>
#include <petscsys.h>

#include "quadrature/tensor_product.h"
#include "basis_on_mesh/basis_on_mesh.h"
#include "spatial_discretization/finite_element_method/01_integrand_stiffness_matrix_laplace.h"

namespace SpatialDiscretization
{

// 1D,2D,3D stiffness matrix of Deformable mesh
template<typename BasisOnMeshType, typename QuadratureType, typename Term>
void AssembleStiffnessMatrix<BasisOnMeshType, QuadratureType, Term>::
setStiffnessMatrix()
{
  const int D = BasisOnMeshType::dim();
  LOG(TRACE) << "setStiffnessMatrix " << D << "D, BasisOnMeshType: " << typeid(BasisOnMeshType).name() << ", QuadratureType: " << typeid(QuadratureType).name();

  // get prefactor value
  const double prefactor = PythonUtility::getOptionDouble(this->specificSettings_, "prefactor", 1.0);

  // define shortcuts for integrator and basis
  typedef Quadrature::TensorProduct<D,QuadratureType> QuadratureDD;
  const int nDofsPerElement = BasisOnMeshType::nDofsPerElement();
  //typedef std::array<std::array<double, nDofsPerElement>, nDofsPerElement> EvaluationsType;
  typedef MathUtility::Matrix<nDofsPerElement,nDofsPerElement> EvaluationsType;
  typedef std::array<
            EvaluationsType,
            QuadratureDD::numberEvaluations()
          > EvaluationsArrayType;     // evaluations[nGP^D][nDofs][nDofs]

  // initialize variables
  std::shared_ptr<PartitionedPetscMat<BasisOnMeshType>> stiffnessMatrix = this->data_.stiffnessMatrix();

  std::shared_ptr<BasisOnMeshType> mesh = std::static_pointer_cast<BasisOnMeshType>(this->data_.mesh());
  mesh->geometryField().startVectorManipulation();

  // initialize values to zero
  int cntr = 1;
  
  LOG(DEBUG) << " nElementsLocal: " << mesh->nElementsLocal();
  
  // loop over elements
  for (element_no_t elementNo = 0; elementNo < mesh->nElementsLocal(); elementNo++)
  {
    std::array<dof_no_t,nDofsPerElement> dofNosLocal = mesh->getElementDofNosLocal(elementNo);

    for (int i=0; i<nDofsPerElement; i++)
    {
      for (int j=0; j<nDofsPerElement; j++)
      {
        VLOG(3) << " initialize stiffnessMatrix entry ( " << dofNosLocal[i] << "," << dofNosLocal[j] << ") (no. " << cntr++ << ")";
        //LOG(DEBUG) << " initialize stiffnessMatrix entry ( " << dofNosLocal[i] << "," << dofNosLocal[j] << ") (no. " << cntr++ << ")";
        stiffnessMatrix->setValue(dofNosLocal[i], dofNosLocal[j], 0, INSERT_VALUES);
      }
    }
  }

  // setup arrays used for integration
  std::array<std::array<double,D>, QuadratureDD::numberEvaluations()> samplingPoints = QuadratureDD::samplingPoints();
  EvaluationsArrayType evaluationsArray{};

  LOG(DEBUG) << "1D integration with " << QuadratureType::numberEvaluations() << " evaluations";
  LOG(DEBUG) << D << "D integration with " << QuadratureDD::numberEvaluations() << " evaluations";
#ifdef DEBUG
  LOG(DEBUG) << "SAMPLING POINTS: ";
  for  (auto value : samplingPoints)
    LOG(DEBUG) << "   " << value;
#endif

  // allow switching between stiffnessMatrix->setValue(... INSERT_VALUES) and ADD_VALUES
  stiffnessMatrix->assembly(MAT_FLUSH_ASSEMBLY);
  
  // fill entries in stiffness matrix
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
      // evaluate function to integrate at samplingPoint
      std::array<double,D> xi = samplingPoints[samplingPointIndex];

      // compute the 3xD jacobian of the parameter space to world space mapping
      auto jacobian = BasisOnMeshType::computeJacobian(geometry, xi);

      VLOG(2) << "samplingPointIndex=" <<samplingPointIndex<< ", xi=" <<xi<< ", geometry: " <<geometry<< ", jac: " <<jacobian;

      // get evaluations of integrand at xi for all (i,j)-dof pairs, integrand is defined in another class
      evaluationsArray[samplingPointIndex]
        = IntegrandStiffnessMatrix<D,EvaluationsType,BasisOnMeshType,Term>::evaluateIntegrand(this->data_,jacobian,xi);

    }  // function evaluations

    // integrate all values for the (i,j) dof pairs at once
    EvaluationsType integratedValues = QuadratureDD::computeIntegral(evaluationsArray);

    // perform integration and add to entry of stiffness matrix
    for (int i=0; i<nDofsPerElement; i++)
    {
      for (int j=0; j<nDofsPerElement; j++)
      {
        // integrate value and set entry in stiffness matrix
        double integratedValue = integratedValues(i,j);
        double value = -prefactor * integratedValue;

        VLOG(2) << "  dof pair (" << i<< "," <<j<< ") dofs (" <<dofNosLocal[i]<< "," <<dofNosLocal[j]<< "), prefactor: " << prefactor << ", integrated value: " <<integratedValue;

        stiffnessMatrix->setValue(dofNosLocal[i], dofNosLocal[j], value, ADD_VALUES);
      }  // j
    }  // i
  }  // elementNo

  // MatAssemblyBegin(stiffnessMatrix,MAT_FINAL_ASSEMBLY) is called after this via data->finalAssembly
}

};
