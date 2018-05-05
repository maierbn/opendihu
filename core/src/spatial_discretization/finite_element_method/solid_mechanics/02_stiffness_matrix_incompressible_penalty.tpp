#include "spatial_discretization/finite_element_method/solid_mechanics/02_stiffness_matrix_incompressible.h"

#include <Python.h>  // has to be the first included header
#include <memory>
#include <petscksp.h>
#include <petscsys.h>
#include <cmath>

#include "semt/Semt.h"
#include "semt/Shortcuts.h"
#include "control/types.h"
#include "spatial_discretization/finite_element_method/solid_mechanics/elasticity_tensor.h"
#include "basis_on_mesh/00_basis_on_mesh_base_dim.h"
#include "utility/python_utility.h"

namespace SpatialDiscretization
{
  
// general implementation for solid mechanics penalty 
template<typename BasisOnMeshType, typename QuadratureType, typename Term>
void FiniteElementMethodStiffnessMatrix<
  BasisOnMeshType,
  QuadratureType,
  Term,
  typename BasisOnMeshType::Mesh,
  Equation::isIncompressible<Term>,
  BasisFunction::isNotMixed<typename BasisOnMeshType::BasisFunction>
>::
setStiffnessMatrix(Mat stiffnessMatrix)
{
  Mat &tangentStiffnessMatrix = (stiffnessMatrix == PETSC_NULL ? this->data_.tangentStiffnessMatrix() : stiffnessMatrix);
  
  // set all non-zero entries
  this->setStiffnessMatrixEntriesForDisplacements(tangentStiffnessMatrix);
  
  // because this is used in nonlinear solver context, assembly has to be called here, not via data->finalAssembly
  PetscErrorCode ierr;
  ierr = MatAssemblyBegin(tangentStiffnessMatrix, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);
  ierr = MatAssemblyEnd(tangentStiffnessMatrix, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);
  
  if (!this->tangentStiffnessMatrixInitialized_)
  {
    VLOG(3) << "tangent stiffness matrix before zero rows columns: " << PetscUtility::getStringMatrix(tangentStiffnessMatrix);
    VLOG(3) << "number dirichletIndices: " << this->dirichletIndices_.size();
    
    // zero rows and columns for which Dirichlet BC is set 
    ierr = MatZeroRowsColumns(tangentStiffnessMatrix, this->dirichletIndices_.size(), this->dirichletIndices_.data(), 1.0, PETSC_IGNORE, PETSC_IGNORE); CHKERRV(ierr);
  
    VLOG(3) << "tangent stiffness matrix after zero rows columns: " << PetscUtility::getStringMatrix(tangentStiffnessMatrix);
    
    // set option that all insert/add operations to new nonzero locations will be discarded. This keeps the nonzero structure forever.
    // (The diagonal entries will be set to different values, but that doesn't matter because the Dirichlet values for updates are 0 and thus independent of the diagonal scaling (d*Δx=0 -> Δx=0 independent of d))
    ierr = MatSetOption(tangentStiffnessMatrix, MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE); CHKERRV(ierr);
    
    this->tangentStiffnessMatrixInitialized_ = true;
  }
}

template<typename BasisOnMeshType, typename QuadratureType, typename Term>
void FiniteElementMethodStiffnessMatrix<
  BasisOnMeshType,
  QuadratureType,
  Term,
  typename BasisOnMeshType::Mesh,
  Equation::isIncompressible<Term>,
  BasisFunction::isNotMixed<typename BasisOnMeshType::BasisFunction>
>::
setFromSolverVariableSolution(Vec &solverSolutionVariable)
{
  // if the computation uses reduced vectors, expand to full vectors
  if (this->data_.computeWithReducedVectors())
  {
    const int D = BasisOnMeshType::dim();
    const int nUnknownsOutputVector = this->data_.mesh()->nDofs() * D;
    
    this->expandVector(solverSolutionVariable, this->data_.displacements().values(), nUnknownsOutputVector);
  }
  else
  {
    VecCopy(solverSolutionVariable, this->data_.displacements().values());
  }
}

template<typename BasisOnMeshType, typename QuadratureType, typename Term>
void FiniteElementMethodStiffnessMatrix<
  BasisOnMeshType,
  QuadratureType,
  Term,
  typename BasisOnMeshType::Mesh,
  Equation::isIncompressible<Term>,
  BasisFunction::isNotMixed<typename BasisOnMeshType::BasisFunction>
>::
evaluateNonlinearFunction(Vec &result)
{
  this->computeInternalMinusExternalVirtualWork(result);
}

template<typename BasisOnMeshType, typename QuadratureType, typename Term>
void FiniteElementMethodStiffnessMatrix<
  BasisOnMeshType,
  QuadratureType,
  Term,
  typename BasisOnMeshType::Mesh,
  Equation::isIncompressible<Term>,
  BasisFunction::isNotMixed<typename BasisOnMeshType::BasisFunction>
>::
preparePressureInterpolation(element_no_t elementNo)
{
  // Do nothing because pressure is not interpolated in penalty formulation. This method is needed to generalize to mixed formulation.
}

template<typename BasisOnMeshType, typename QuadratureType, typename Term>
double FiniteElementMethodStiffnessMatrix<
  BasisOnMeshType,
  QuadratureType,
  Term,
  typename BasisOnMeshType::Mesh,
  Equation::isIncompressible<Term>,
  BasisFunction::isNotMixed<typename BasisOnMeshType::BasisFunction>
>::
getPressure(double deformationGradientDeterminant, VecD<BasisOnMeshType::dim()> xi, double &artificialPressureTilde)
{
  // artifical pressure p
  const double artificialPressure = this->computeArtificialPressure(deformationGradientDeterminant, artificialPressureTilde);
 
  return artificialPressure;
}

template<typename BasisOnMeshType, typename QuadratureType, typename Term>
const int FiniteElementMethodStiffnessMatrix<
  BasisOnMeshType,
  QuadratureType,
  Term,
  typename BasisOnMeshType::Mesh,
  Equation::isIncompressible<Term>,
  BasisFunction::isNotMixed<typename BasisOnMeshType::BasisFunction>
>::
nUnknowns()
{
  const int D = BasisOnMeshType::dim();
  const int nUnknowns = this->data_.mesh()->nDofs() * D;
  return nUnknowns;
}

};    // namespace