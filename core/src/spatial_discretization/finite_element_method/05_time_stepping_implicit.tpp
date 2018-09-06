#include "spatial_discretization/finite_element_method/05_time_stepping.h"

#include <Python.h>
#include <iostream>
#include <petscmat.h>
#include <petscksp.h>
#include <vector>
#include <numeric>

#include "easylogging++.h"

#include "control/types.h"
#include "utility/python_utility.h"
#include "utility/petsc_utility.h"
#include "solver/solver_manager.h"
#include "solver/linear.h"


namespace SpatialDiscretization
{

/* 
template<typename FunctionSpaceType, typename QuadratureType, typename Term>
void FiniteElementMethodTimeStepping<FunctionSpaceType, QuadratureType, Term>::
evaluateTimesteppingRightHandSideImplicit(Vec &input, Vec &output, int timeStepNo, double currentTime)
{
}
*/

template<typename FunctionSpaceType, typename QuadratureType, typename Term>
void FiniteElementMethodTimeStepping<FunctionSpaceType, QuadratureType, Term>::
setInverseLumpedMassMatrix()
{
  LOG(TRACE) << "setInverseLumpedMassMatrix";

  std::shared_ptr<FunctionSpaceType> functionSpace = std::static_pointer_cast<FunctionSpaceType>(this->data_.functionSpace());
  Mat &inverseLumpedMassMatrix = this->data_.inverseLumpedMassMatrix()->valuesGlobal();
  Mat &massMatrix = this->data_.massMatrix()->valuesGlobal();
  //this->data_.massMatrix()->assembly(MAT_FINAL_ASSEMBLY);
  
  PetscErrorCode ierr;
  
  PetscInt nRows, nColumns;
  ierr = MatGetSize(massMatrix,&nRows,&nColumns); CHKERRV(ierr);
  VLOG(1) << "massMatrix nRows " << nRows << " nColumns " << nColumns;
     
  std::shared_ptr<PartitionedPetscVec<FunctionSpaceType,1>> rowSum = std::make_shared<PartitionedPetscVec<FunctionSpaceType,1>>(functionSpace->meshPartition(), "rowSum");

  // In case of linear and bilinear basis functions
  // store the sum of each row of the matrix in the vector rowSum
  ierr = MatGetRowSum(massMatrix, rowSum->valuesGlobal()); CHKERRV(ierr);

  // for the inverse matrix, replace each entry in rowSum by its reciprocal
  ierr = VecReciprocal(rowSum->valuesGlobal()); CHKERRV(ierr);

  // set the values on the diagonal
  ierr = MatDiagonalSet(inverseLumpedMassMatrix, rowSum->valuesGlobal(), INSERT_VALUES); CHKERRV(ierr);

  this->data_.inverseLumpedMassMatrix()->assembly(MAT_FINAL_ASSEMBLY);

  VLOG(2) << *this->data_.inverseLumpedMassMatrix();
}

} // namespace SpatialDiscretization
