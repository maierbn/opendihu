#include "spatial_discretization/finite_element_method/01_matrix.h"

#include <Python.h>  // has to be the first included header
#include <memory>
#include <vector>
#include <petscsys.h>
#include <array>

#include "quadrature/tensor_product.h"
#include "function_space/function_space.h"
#include "spatial_discretization/finite_element_method/integrand/integrand_stiffness_matrix_laplace.h"
#include "spatial_discretization/finite_element_method/integrand/integrand_stiffness_matrix_linear_elasticity.h"
#include "control/types.h"


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
  typedef MathUtility::Matrix<nUnknownsPerElement,nUnknownsPerElement,double_v_t> EvaluationsType;
  typedef std::array<
            EvaluationsType,
            QuadratureDD::numberEvaluations()
          > EvaluationsArrayType;     // evaluations[nGP^D][nDofs][nDofs]

  // setup arrays used for integration
  std::array<std::array<double,D>, QuadratureDD::numberEvaluations()> samplingPoints = QuadratureDD::samplingPoints();
  EvaluationsArrayType evaluationsArray{};

  LOG(DEBUG) << "1D integration with " << QuadratureType::numberEvaluations() << " evaluations";
  LOG(DEBUG) << D << "D integration with " << QuadratureDD::numberEvaluations() << " evaluations";

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

  const element_no_t nElementsLocal = functionSpace->nElementsLocal();
  LOG(DEBUG) << " nElementsLocal: " << nElementsLocal;

  // initialize values to zero
  // loop over elements, always 4 elements at once using the vectorized functions
  for (int elementNoLocal = 0; elementNoLocal < nElementsLocal; elementNoLocal += nVcComponents)
  {

#ifdef USE_VC
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

    std::array<dof_no_v_t,nDofsPerElement> dofNosLocal = functionSpace->getElementDofNosLocal(elementNoLocalv);

    for (int i = 0; i < nDofsPerElement; i++)
    {
      for (int j = 0; j < nDofsPerElement; j++)
      {
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

  // allow switching between stiffnessMatrix->setValue(... INSERT_VALUES) and ADD_VALUES
  stiffnessMatrix->assembly(MAT_FLUSH_ASSEMBLY);
  
  double progress = 0;

  // fill entries in stiffness matrix
  // loop over elements, always 4 elements at once using the vectorized functions
  for (int elementNoLocal = 0; elementNoLocal < nElementsLocal; elementNoLocal += nVcComponents)
  {

#ifdef USE_VC
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
    std::array<dof_no_v_t,nDofsPerElement> dofNosLocal = functionSpace->getElementDofNosLocal(elementNoLocalv);

    VLOG(2) << "element " << elementNoLocalv;

    // get geometry field (which are the node positions for Lagrange basis and node positions and derivatives for Hermite)
    std::array<Vec3_v_t,FunctionSpaceType::nDofsPerElement()> geometry;
    functionSpace->getElementGeometry(elementNoLocalv, geometry);

    // compute integral
    for (unsigned int samplingPointIndex = 0; samplingPointIndex < samplingPoints.size(); samplingPointIndex++)
    {
      // evaluate function to integrate at samplingPoint
      std::array<double,D> xi = samplingPoints[samplingPointIndex];

      // compute the 3xD jacobian of the parameter space to world space mapping
      std::array<Vec3_v_t,D> jacobian = FunctionSpaceType::computeJacobian(geometry, xi);

      VLOG(2) << "samplingPointIndex=" << samplingPointIndex<< ", xi=" <<xi<< ", geometry: " <<geometry<< ", jac: " <<jacobian;

      const double_v_t prefactor = this->prefactor_.value(elementNoLocalv);

      // get evaluations of integrand at xi for all (i,j)-dof pairs, integrand is defined in another class
      // gradPhi[j](xi)^T * T * gradPhi[k](xi)
      evaluationsArray[samplingPointIndex]
        = prefactor * IntegrandStiffnessMatrix<D,EvaluationsType,FunctionSpaceType,nComponents,double_v_t,dof_no_v_t,Term>::
          evaluateIntegrand(this->data_, jacobian, elementNoLocalv, xi);

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
            double_v_t integratedValue = integratedValues(i*nComponents + rowComponentNo, j*nComponents + columnComponentNo);
            double_v_t value = -integratedValue;
            int componentNo = rowComponentNo*nComponents + columnComponentNo;

            VLOG(2) << "  dof pair (" << i<< "," <<j<< ") dofs (" << dofNosLocal[i]<< "," << dofNosLocal[j]<< "), "
              << "component (" << rowComponentNo << "," << columnComponentNo << "), " << componentNo
              << ", integrated value: " <<integratedValue;

            // get local dof no
            dof_no_v_t dofINoLocal = dofNosLocal[i];
            dof_no_v_t dofJNoLocal = dofNosLocal[j];

            // add the entry in the stiffness matrix, for all dofs of the vectorized values at once,
            // i.e. K_dofINoLocal[0],dofJNoLocal[0] = value[0]
            // i.e. K_dofINoLocal[1],dofJNoLocal[1] = value[1], etc.
            // Note that K_dofINoLocal[0],dofJNoLocal[1] would be potentially zero, the contributions are considered element-wise
            stiffnessMatrix->setValue(componentNo, dofINoLocal, dofJNoLocal, value, ADD_VALUES);
          }
        }
      }  // j
    }  // i
  }  // elementNoLocalv

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
