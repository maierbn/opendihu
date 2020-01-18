#include "spatial_discretization/finite_element_method/02_boundary_conditions.h"

#include <Python.h>
#include <memory>
#include <vector>
#include <petscsys.h>

#include "quadrature/tensor_product.h"
#include "utility/vector_operators.h"

namespace SpatialDiscretization
{

template<typename FunctionSpaceType,typename QuadratureType,int nComponents,typename Term,typename Dummy>
void BoundaryConditions<FunctionSpaceType,QuadratureType,nComponents,Term,Dummy>::
setBoundaryConditionHandlingEnabled(bool boundaryConditionHandlingEnabled)
{
  boundaryConditionHandlingEnabled_ = boundaryConditionHandlingEnabled;
}

template<typename FunctionSpaceType,typename QuadratureType,int nComponents,typename Term,typename Dummy>
void BoundaryConditions<FunctionSpaceType,QuadratureType,nComponents,Term,Dummy>::
setDirichletBoundaryConditions(std::shared_ptr<DirichletBoundaryConditions<FunctionSpaceType,nComponents>> dirichletBoundaryConditions)
{
  this->dirichletBoundaryConditions_ = dirichletBoundaryConditions;
}

template<typename FunctionSpaceType,typename QuadratureType,int nComponents,typename Term,typename Dummy>
void BoundaryConditions<FunctionSpaceType,QuadratureType,nComponents,Term,Dummy>::
setNeumannBoundaryConditions(std::shared_ptr<NeumannBoundaryConditions<FunctionSpaceType,QuadratureType,nComponents>> neumannBoundaryConditions)
{
  LOG(DEBUG) << "set Neumann boundary conditions";
  this->neumannBoundaryConditions_ = neumannBoundaryConditions;
  neumannBoundaryConditionsApplied_ = false;
}

template<typename FunctionSpaceType,typename QuadratureType,int nComponents,typename Term,typename Dummy>
void BoundaryConditions<FunctionSpaceType,QuadratureType,nComponents,Term,Dummy>::
reset()
{
  LOG(DEBUG) << "delete dirichlet boundary conditions object";
  this->dirichletBoundaryConditions_ = nullptr;
  this->systemMatrixAlreadySet_ = false;
}

template<typename FunctionSpaceType,typename QuadratureType,int nComponents,typename Term,typename Dummy>
void BoundaryConditions<FunctionSpaceType,QuadratureType,nComponents,Term,Dummy>::
applyBoundaryConditions()
{
  if (!boundaryConditionHandlingEnabled_)
  {
    if (this->specificSettings_.hasKey("dirichletBoundaryConditions"))
    {
      LOG(WARNING) << "You have specified dirichlet boundary conditions in FiniteElementMethod via the key \"dirichletBoundaryConditions\". "
        << "They are not used here, e.g. because the FiniteElementMethod is wrapped by a time stepping scheme. Consider only setting dirichlet boundary conditions in the time stepping scheme.";
    }

    VLOG(1) << "do not handle boundary conditions in finite element method, because boundaryConditionHandlingEnabled=false";
    return;
  }

  LOG(TRACE) << "FiniteElementMethod::applyBoundaryConditions";

  if (VLOG_IS_ON(4))
  {
    VLOG(4) << "Finite element data before applyBoundaryConditions";
    this->data_.print();
  }

  // handle Neumann boundary conditions
  applyNeumannBoundaryConditions();

  // handle Dirichlet boundary conditions
  applyDirichletBoundaryConditions();

  if (VLOG_IS_ON(4))
  {
    VLOG(4) << "Finite element data after applyBoundaryConditions";
    this->data_.print();
  }
}

template<typename FunctionSpaceType,typename QuadratureType,int nComponents,typename Term,typename Dummy>
void BoundaryConditions<FunctionSpaceType,QuadratureType,nComponents,Term,Dummy>::
applyNeumannBoundaryConditions()
{
  // handle Neumann boundary conditions, only change the rhs once, not in every timestep, if this solver is in a timestepping scheme
  if (neumannBoundaryConditions_ == nullptr)
  {
    LOG(DEBUG) << "no Neumann boundary conditions are present, create object";
    neumannBoundaryConditions_ = std::make_shared<NeumannBoundaryConditions<FunctionSpaceType,QuadratureType,nComponents>>(this->context_);
    neumannBoundaryConditions_->initialize(this->specificSettings_, this->data_.functionSpace(), "neumannBoundaryConditions");
    this->data_.setNegativeRightHandSideNeumannBoundaryConditions(neumannBoundaryConditions_->rhs());
  }

  if (!neumannBoundaryConditionsApplied_)
  {
    neumannBoundaryConditionsApplied_ = true;

    LOG(DEBUG) << "neumann BC rhs: " << *neumannBoundaryConditions_->rhs();
    LOG(DEBUG) << "rhs: " << *this->data_.rightHandSide();

    // add rhs, rightHandSide += -1 * rhs
    PetscErrorCode ierr;
    ierr = VecAXPY(this->data_.rightHandSide()->valuesGlobal(), -1, neumannBoundaryConditions_->rhs()->valuesGlobal()); CHKERRV(ierr);
  }
}

template<typename FunctionSpaceType,typename QuadratureType,int nComponents,typename Term,typename Dummy>
void BoundaryConditions<FunctionSpaceType,QuadratureType,nComponents,Term,Dummy>::
applyDirichletBoundaryConditions()
{
  // handle Dirichlet boundary conditions
  if (dirichletBoundaryConditions_ == nullptr)
  {
    LOG(DEBUG) << "no Dirichlet boundary conditions are present, create object";
    dirichletBoundaryConditions_ = std::make_shared<DirichletBoundaryConditions<FunctionSpaceType,nComponents>>(this->context_);
    dirichletBoundaryConditions_->initialize(this->specificSettings_, this->data_.functionSpace(), "dirichletBoundaryConditions");
  }
  else
  {
    LOG(DEBUG) << "dirichlet BC object already exists";
  }

  // if the option to use the bc values from solution is set
  if (this->updatePrescribedValuesFromSolution_)
  {
    LOG(DEBUG) << "updatePrescribedValuesFromSolution";

    // update the prescribed boundary condition values
    dirichletBoundaryConditions_->updatePrescribedValuesFromSolution(this->data_.solution());

    // reset right hand side
    this->setRightHandSide();
  }

  // get abbreviations
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nComponents>> rightHandSide = this->data_.rightHandSide();
  std::shared_ptr<PartitionedPetscMat<FunctionSpaceType>> stiffnessMatrix = this->data_.stiffnessMatrix();
  std::shared_ptr<PartitionedPetscMat<FunctionSpaceType>> stiffnessMatrixWithoutBc = this->data_.stiffnessMatrixWithoutBc();

  // apply the boundary conditions in stiffness matrix
  // (set bc rows and columns of stiffnessMatrix to 0 and diagonal to 1), also add terms with matrix entries to rhs, for the reading of matrix entries, stiffnessMatrixWithoutBc is used.
  LOG(DEBUG) << "call applyInSystemMatrix from applyBoundaryConditions, this->systemMatrixAlreadySet: " << this->systemMatrixAlreadySet_;
  dirichletBoundaryConditions_->applyInSystemMatrix(stiffnessMatrixWithoutBc, stiffnessMatrix, rightHandSide, this->systemMatrixAlreadySet_);
  this->systemMatrixAlreadySet_ = true;

  // set prescribed values in rhs
  dirichletBoundaryConditions_->applyInRightHandSide(rightHandSide, rightHandSide);
}

} // namespace
