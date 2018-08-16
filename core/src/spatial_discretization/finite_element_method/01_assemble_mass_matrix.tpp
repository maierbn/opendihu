#include "spatial_discretization/finite_element_method/01_assemble_finite_element_matrix.h"

#include <Python.h>
#include <memory>
#include <petscksp.h>
#include <petscsys.h>

#include "mesh/structured_regular_fixed.h"
#include "basis_function/lagrange.h"

#include <Python.h>
#include <memory>

#include "spatial_discretization/spatial_discretization.h"
#include "spatial_discretization/finite_element_method/03_integrand_rhs.h"
//#include "time_stepping_scheme/discretizable_in_time.h"
#include "control/runnable.h"
#include "control/dihu_context.h"
#include "utility/petsc_utility.h"
#include "data_management/finite_elements.h"
#include "equation/laplace.h"
#include "equation/poisson.h"
#include "equation/type_traits.h"
#include "mesh/mesh.h"
#include "control/types.h"


namespace SpatialDiscretization
{

// 1D,2D,3D mass matrix of Deformable mesh
template<typename BasisOnMeshType, typename QuadratureType, typename Term>
void AssembleFiniteElementMatrix<BasisOnMeshType, QuadratureType, Term>::
setMassMatrix()
{
  // check if matrix discretization matrix exists
  if (!this->data_.massMatrixInitialized())
  {
    this->data_.initializeMassMatrix();

    const int D = BasisOnMeshType::dim();
    LOG(INFO)<<"setMassMatrix " << D << "D";

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
    PetscErrorCode ierr;
    Mat &massMatrix = this->data_.massMatrix();

    std::shared_ptr<BasisOnMeshType> mesh = std::static_pointer_cast<BasisOnMeshType>(this->data_.mesh());

    // initialize values to zero
    // loop over elements
    for (element_no_t elementNo = 0; elementNo < mesh->nElements(); elementNo++)
    {
      auto dof = mesh->getElementDofNos(elementNo);

      for (int i=0; i<nDofsPerElement; i++)
      {
        for (int j=0; j<nDofsPerElement; j++)
        {
          ierr = MatSetValue(massMatrix, dof[i], dof[j], 0, INSERT_VALUES); CHKERRV(ierr);
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

    // set entries in massMatrix
    // loop over elements
    for (element_no_t elementNo = 0; elementNo < mesh->nElements(); elementNo++)
    {
      // get indices of element-local dofs
      auto dof = mesh->getElementDofNos(elementNo);

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

          ierr = MatSetValue(massMatrix, dof[i], dof[j], integratedValue, ADD_VALUES); CHKERRV(ierr);
        }  // j
      }  // i
    }  // elementNo

    ierr = MatAssemblyBegin(massMatrix, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);
    ierr = MatAssemblyEnd(massMatrix, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);
    
    /*
    dof_no_t n = mesh->nDofs();
    int nEntries=n;       
    PetscScalar v=0.0;  
    for (PetscInt i=0;i<nEntries;i++)
      for(PetscInt j=0;j<nEntries;j++)
      {
        ierr=MatGetValues(massMatrix,1,&i,1,&j,&v);
        LOG(DEBUG)<<"val_get: massMatrix assembly"<< v; 
      }
      */
  }
}


};    // namespace