#include "spatial_discretization/finite_element_method/01_matrix.h"

#include <Python.h>  // has to be the first included header
#include <memory>
#include <vector>
#include <petscsys.h>
#include <Vc/Vc>

#include "quadrature/tensor_product.h"
#include "function_space/function_space.h"
#include "spatial_discretization/finite_element_method/integrand/integrand_stiffness_matrix_laplace.h"
#include "spatial_discretization/finite_element_method/integrand/integrand_stiffness_matrix_linear_elasticity.h"

namespace SpatialDiscretization
{
// 1D,2D,3D stiffness matrix of Deformable mesh
template<typename FunctionSpaceType,typename QuadratureType,int nComponents,typename Term,typename Dummy1, typename Dummy2, typename Dummy3>
void FiniteElementMethodMatrix<FunctionSpaceType,QuadratureType,nComponents,Term,Dummy1,Dummy2,Dummy3>::
setStiffnessMatrix()
{
  const int D = FunctionSpaceType::dim();
  LOG(TRACE) << "setStiffnessMatrix " << D << "D using integration, FunctionSpaceType: " << StringUtility::demangle(typeid(FunctionSpaceType).name())
    << ", QuadratureType: " << StringUtility::demangle(typeid(QuadratureType).name());

  // define shortcuts for integrator and basis
  typedef Quadrature::TensorProduct<D,QuadratureType> QuadratureDD;
  const int nDofsPerElement = FunctionSpaceType::nDofsPerElement();
  const int nUnknownsPerElement = nDofsPerElement*nComponents;
  typedef MathUtility::Matrix<nUnknownsPerElement,nUnknownsPerElement> EvaluationsType;
  typedef std::array<
            EvaluationsType,
            QuadratureDD::numberEvaluations()
          > EvaluationsArrayType;     // evaluations[nGP^D][nDofs][nDofs]

  // initialize variables
  std::shared_ptr<PartitionedPetscMat<FunctionSpaceType>> stiffnessMatrix = this->data_.stiffnessMatrix();

  std::shared_ptr<FunctionSpaceType> functionSpace = std::static_pointer_cast<FunctionSpaceType>(this->data_.functionSpace());
  functionSpace->geometryField().setRepresentationGlobal();
  functionSpace->geometryField().startGhostManipulation();   // ensure that local ghost values of geometry field are set

  bool outputAssemble3DStiffnessMatrixHere = false;
  if (outputAssemble3DStiffnessMatrix_)
  {
    LOG(INFO) << "Compute stiffness matrix for " << D << "D problem with " << functionSpace->nDofsGlobal() << " global dofs.";
    outputAssemble3DStiffnessMatrix_ = false;
    outputAssemble3DStiffnessMatrixHere = true;
  }

  // initialize values to zero
  //int cntr = 1;

  LOG(DEBUG) << " nElementsLocal: " << functionSpace->nElementsLocal();

  // loop over elements
  for (element_no_t elementNoLocal = 0; elementNoLocal < functionSpace->nElementsLocal(); elementNoLocal++)
  {
    std::array<dof_no_t,nDofsPerElement> dofNosLocal = functionSpace->getElementDofNosLocal(elementNoLocal);
/*
    bool show = false;
    for (int i = 0; i < nDofsPerElement; i++)
    {
      if (dofNosLocal[i] == 27)
        show = true;
    }

    if (show)
    {
      for (int i = 0; i < nDofsPerElement; i++)
      {
        Vec3 nodePosition0 = functionSpace->getGeometry(dofNosLocal[i]);

        VLOG(1) << "node_position_dof" << dofNosLocal[i] << " = " << nodePosition0;
      }
    }
*/
    for (int i = 0; i < nDofsPerElement; i++)
    {
      for (int j = 0; j < nDofsPerElement; j++)
      {
/*
        if (show)
        {
          VLOG(1) << " initialize stiffnessMatrix entry for element " << elementNoLocal << " elementalDofs (" << i << "," << j << "), "
            << "localDofs " << dofNosLocal[i] << "," << dofNosLocal[j] << ") (entry no. " << cntr++ << ")";
        }
*/
        // loop over components (1,...,D for solid mechanics)
        for (int rowComponentNo = 0; rowComponentNo < nComponents; rowComponentNo++)
        {
          for (int columnComponentNo = 0; columnComponentNo < nComponents; columnComponentNo++)
          {
            int componentNo = rowComponentNo*nComponents + columnComponentNo;

            //LOG(DEBUG) << " initialize stiffnessMatrix entry ( " << dofNosLocal[i] << "," << dofNosLocal[j] << ") (no. " << cntr++ << ")";
            stiffnessMatrix->setValue(componentNo, dofNosLocal[i], dofNosLocal[j], 0, INSERT_VALUES);
          }
        }
      }
    }
  }

  // setup arrays used for integration
  Vc::array<std::array<double,D>, QuadratureDD::numberEvaluations()> samplingPoints = QuadratureDD::samplingPoints();
  EvaluationsArrayType evaluationsArray{};

  LOG(DEBUG) << "1D integration with " << QuadratureType::numberEvaluations() << " evaluations";
  LOG(DEBUG) << D << "D integration with " << QuadratureDD::numberEvaluations() << " evaluations";

  // allow switching between stiffnessMatrix->setValue(... INSERT_VALUES) and ADD_VALUES
  stiffnessMatrix->assembly(MAT_FLUSH_ASSEMBLY);
  
  double progress = 0;
  element_no_t nElementsLocal = functionSpace->nElementsLocal();

  // fill entries in stiffness matrix
  // loop over elements
  for (element_no_t elementNoLocal = 0; elementNoLocal < nElementsLocal; elementNoLocal++)
  {
    if (outputAssemble3DStiffnessMatrixHere && this->context_.ownRankNoCommWorld() == 0)
    {
      double newProgress = (double)elementNoLocal / nElementsLocal;
      if (int(newProgress*10) != int(progress*10))
      {
        std::cout << "\b\b\b\b" << int(newProgress*100) << "%" << std::flush;
      }
      progress = newProgress;
    }

    // get indices of element-local dofs
    std::array<dof_no_t,nDofsPerElement> dofNosLocal = functionSpace->getElementDofNosLocal(elementNoLocal);

    VLOG(2) << "element " << elementNoLocal;

    // get geometry field (which are the node positions for Lagrange basis and node positions and derivatives for Hermite)
    std::array<Vec3,FunctionSpaceType::nDofsPerElement()> geometry;
    functionSpace->getElementGeometry(elementNoLocal, geometry);

    // compute integral
    for (unsigned int samplingPointIndex = 0; samplingPointIndex < samplingPoints.size(); samplingPointIndex++)
    {
      // evaluate function to integrate at samplingPoint
      std::array<double,D> xi = samplingPoints[samplingPointIndex];

      // compute the 3xD jacobian of the parameter space to world space mapping
      std::array<Vec3,D> jacobian = FunctionSpaceType::computeJacobian(geometry, xi);

      VLOG(2) << "samplingPointIndex=" << samplingPointIndex<< ", xi=" <<xi<< ", geometry: " <<geometry<< ", jac: " <<jacobian;

      const double prefactor = this->prefactor_.value(elementNoLocal);

      // get evaluations of integrand at xi for all (i,j)-dof pairs, integrand is defined in another class
      // gradPhi[j](xi)^T * T * gradPhi[k](xi)
      evaluationsArray[samplingPointIndex]
        = prefactor * IntegrandStiffnessMatrix<D,EvaluationsType,FunctionSpaceType,nComponents,Term>::
          evaluateIntegrand(this->data_, jacobian, elementNoLocal, xi);

          /*
      for (int i = 0; i < nDofsPerElement; i++)
      {
        for (int j = 0; j < nDofsPerElement; j++)
        {
          // loop over components (1,...,D for solid mechanics)
          for (int rowComponentNo = 0; rowComponentNo < nComponents; rowComponentNo++)
          {
            for (int columnComponentNo = 0; columnComponentNo < nComponents; columnComponentNo++)
            {
              // integrate value and set entry in stiffness matrix
              double evaluatedValue = evaluationsArray[samplingPointIndex](i*nComponents + rowComponentNo, j*nComponents + columnComponentNo);

              if (rowComponentNo != columnComponentNo)
                LOG(DEBUG) << rowComponentNo << columnComponentNo << ": evaluatedValue: " << evaluatedValue
                  << " at (" << i*nComponents + rowComponentNo << "," << j*nComponents + columnComponentNo << ")";
            }
          }
        }
      }
      */

    }  // function evaluations

    // integrate all values for the (i,j) dof pairs at once
    EvaluationsType integratedValues = QuadratureDD::computeIntegral(evaluationsArray);

    // perform integration and add to entry of stiffness matrix
    for (int i = 0; i < nDofsPerElement; i++)
    {
      for (int j = 0; j < nDofsPerElement; j++)
      {
        // loop over components (1,...,D for solid mechanics)
        for (int rowComponentNo = 0; rowComponentNo < nComponents; rowComponentNo++)
        {
          for (int columnComponentNo = 0; columnComponentNo < nComponents; columnComponentNo++)
          {
            // integrate value and set entry in stiffness matrix
            double integratedValue = integratedValues(i*nComponents + rowComponentNo, j*nComponents + columnComponentNo);
            double value = -integratedValue;
            int componentNo = rowComponentNo*nComponents + columnComponentNo;

            VLOG(2) << "  dof pair (" << i<< "," <<j<< ") dofs (" << dofNosLocal[i]<< "," << dofNosLocal[j]<< "), "
              << "component (" << rowComponentNo << "," << columnComponentNo << "), " << componentNo
              << ", integrated value: " <<integratedValue;

            stiffnessMatrix->setValue(componentNo, dofNosLocal[i], dofNosLocal[j], value, ADD_VALUES);
          }
        }
      }  // j
    }  // i
  }  // elementNoLocal

  if (outputAssemble3DStiffnessMatrixHere && this->context_.ownRankNoCommWorld() == 0)
  {
    std::cout << "\b\b\b\bparallel assembly..." << std::flush;
  }

  stiffnessMatrix->assembly(MAT_FINAL_ASSEMBLY);

  if (outputAssemble3DStiffnessMatrixHere && this->context_.ownRankNoCommWorld() == 0)
  {
    std::cout << std::string(100,'\b') << "done.                       " << std::endl;
  }

#ifndef NDEBUG
  //LOG(DEBUG) << "stiffnessMatrix:";
  //LOG(DEBUG) << *stiffnessMatrix;
#endif

}

}  // namespace
