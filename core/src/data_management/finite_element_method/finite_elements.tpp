#include "data_management/finite_element_method/finite_elements.h"

#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>
#include <memory>

#include "easylogging++.h"

#include "utility/python_utility.h"
#include "control/dihu_context.h"
#include "utility/petsc_utility.h"
#include "function_space/function_space.h"
#include "mesh/unstructured_deformable.h"
#include "basis_function/hermite.h"
#include "partition/partitioned_petsc_mat/partitioned_petsc_mat.h"

namespace Data
{

//! constructor
template<typename FunctionSpaceType,int nComponents>
FiniteElements<FunctionSpaceType,nComponents,Equation::Dynamic::AnisotropicDiffusion>::
FiniteElements(DihuContext context) :
  FiniteElementsBase<FunctionSpaceType,nComponents>(context),
  DiffusionTensorConstant<FunctionSpaceType>(context.getPythonConfig())
{
}

//! constructor
template<typename FunctionSpaceType,int nComponents>
FiniteElements<FunctionSpaceType,nComponents,Equation::Dynamic::DirectionalDiffusion>::
FiniteElements(DihuContext context) :
  FiniteElementsBase<FunctionSpaceType,nComponents>(context),
  DiffusionTensorDirectional<FunctionSpaceType>(context.getPythonConfig())
{
}

template<typename FunctionSpaceType,int nComponents>
void FiniteElements<FunctionSpaceType,nComponents,Equation::Dynamic::AnisotropicDiffusion>::
initialize()
{
  LOG(DEBUG) << "Data::FiniteElements::initialize";

  FiniteElementsBase<FunctionSpaceType,nComponents>::initialize();

  // set up diffusion tensor if there is any
  DiffusionTensorConstant<FunctionSpaceType>::initialize();
}

template<typename FunctionSpaceType,int nComponents>
void FiniteElements<FunctionSpaceType,nComponents,Equation::Dynamic::DirectionalDiffusion>::
initialize(std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>> direction,
           std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>> spatiallyVaryingPrefactor,
           bool useAdditionalDiffusionTensor)
{
  FiniteElementsBase<FunctionSpaceType,nComponents>::initialize();

  // set up diffusion tensor, initialize with given direction field
  DiffusionTensorDirectional<FunctionSpaceType>::initialize(direction, spatiallyVaryingPrefactor, useAdditionalDiffusionTensor);
}

//! initialize, store the reference geometry as copy of the current geometry
template<typename FunctionSpaceType,int nComponents,typename Term>
void FiniteElements<FunctionSpaceType,nComponents,Term,Equation::isLinearElasticity<Term>>::
initialize()
{
  LinearStiffness<FunctionSpaceType,nComponents>::initialize();

  referenceGeometry_ = this->functionSpace_->template createFieldVariable<3>("referenceGeometry");
  referenceGeometry_->setValues(this->functionSpace_->geometryField());
}

//! update the geometry of the mesh and function space with the displacements
template<typename FunctionSpaceType,int nComponents,typename Term>
void FiniteElements<FunctionSpaceType,nComponents,Term,Equation::isLinearElasticity<Term>>::
updateGeometry()
{
  PetscErrorCode ierr;

  // w = alpha * x + y, VecWAXPY(w, alpha, x, y)
  ierr = VecWAXPY(this->functionSpace_->geometryField().valuesGlobal(), 1, this->referenceGeometry_->valuesGlobal(), this->solution()->valuesGlobal()); CHKERRV(ierr);
}

//! compute the linear strain field epsilon
template<typename FunctionSpaceType,int nComponents,typename Term>
void FiniteElements<FunctionSpaceType,nComponents,Term,Equation::isLinearElasticity<Term>>::
computeStrain(std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nComponents*nComponents>> strain)
{
  // compute strain as epsilon_ab = 1/2(u_Lb dphi_L(x)/dx_a + u_Ma dphi_M(x)/dx_b)
  // compute strain as epsilon_ab = 1/2*u_Lb dphi_L(x)/dx_a + 1/2*u_La dphi_L(x)/dx_b)

  const int D = FunctionSpaceType::dim();

  const int nDofsPerElement = FunctionSpaceType::nDofsPerElement();

  // initialize variables
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nComponents>> displacements = this->solution();

  // We set values at nodes, without any integration. The gradient dphi_L/dx_a is potentially discontinous, therefore we loop over elements.
  // The value at each node is set multiple times, once from each adjacent element, the last value wins.
  // It is important to not communicate ghost values here, because we do not add values.

  // loop over local elements
  for (element_no_t elementNoLocal = 0; elementNoLocal < this->functionSpace_->nElementsLocal(); elementNoLocal++)
  {
    // get indices of element-local dofs
    std::array<dof_no_t,nDofsPerElement> dofNosLocal = this->functionSpace_->getElementDofNosLocal(elementNoLocal);

    // get geometry field (which are the node positions for Lagrange basis and node positions and derivatives for Hermite)
    std::array<Vec3,nDofsPerElement> geometryValues;
    this->functionSpace_->getElementGeometry(elementNoLocal, geometryValues);

    std::array<VecD<nComponents>,nDofsPerElement> displacementValues;
    displacements->getElementValues(elementNoLocal, displacementValues);

    // loop over dofs in element, where to compute the gradient
    for (int dofIndexL = 0; dofIndexL < nDofsPerElement; dofIndexL++)
    {
      dof_no_t dofNoLocal = dofNosLocal[dofIndexL];
      std::array<double,D> xi{};

      // set xi to dofIndexL
      for (int i = 0; i < D; i++)
      {
        if (i == 0)
          xi[i] = double(dofIndexL % 2);
        else if (i == 1)
          xi[i] = double((dofIndexL % 4) / 2);
        else if (i == 2)
          xi[i] = double(dofIndexL / 4);
      }

      // compute gradient values in parameter space
      std::array<VecD<D>,nDofsPerElement> gradPhi = this->functionSpace_->getGradPhi(xi);

      // compute inverse jacobian
      Tensor2<D> inverseJacobian = this->functionSpace_->getInverseJacobian(geometryValues, elementNoLocal, xi);

      VecD<D> dPhiL_dX = inverseJacobian * gradPhi[dofIndexL];

      std::array<double,nComponents*nComponents> epsilon{};
      for (int a = 0; a < D; a++)
      {
        for (int b = 0; b < D; b++)
        {
          epsilon[a*D + b] = 0.5 * (displacementValues[dofIndexL][b] * dPhiL_dX[a] + displacementValues[dofIndexL][a] * dPhiL_dX[b]);
        }
      }

      // set value at node
      strain->setValue(dofNoLocal, epsilon, INSERT_VALUES);

      LOG(DEBUG) << "compute strain dof " << dofIndexL << ", displacements: " << displacementValues[dofIndexL]
        << ", strain: " << epsilon;

    } // dofIndexL
  }  // elementNoLocal
}

} // namespace Data
