#include "finite_element_method_time_stepping/05_finite_element_method_time_stepping.h"

#include <Python.h>
#include <iostream>
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

template<typename BasisOnMeshType, typename QuadratureType, typename Term>
FiniteElementMethodTimeStepping<BasisOnMeshType, QuadratureType, Term>::
FiniteElementMethodTimeStepping(DihuContext context)
  : FiniteElementMethodBaseRhs<BasisOnMeshType, QuadratureType, Term>(context),
  DiscretizableInTime(SolutionVectorMapping(true))
{
  // the solutionVectorMapping_ object stores the information which range of values of the solution will be further used
  // in methods that use the result of this method, e.g. in operator splittings. Since there are no internal values
  // in this FEM, set the range to all values.
  solutionVectorMapping_.setOutputRange(0, this->data_.mesh()->nNodes());
}

template<typename BasisOnMeshType, typename QuadratureType, typename Term>
void FiniteElementMethodTimeStepping<BasisOnMeshType, QuadratureType, Term>::
initialize()
{
  LOG(DEBUG) << "FiniteElementMethodTimeStepping::initialize";
  
  FiniteElementMethodBase<BasisOnMeshType,QuadratureType,Term>::initialize();
}

template<typename BasisOnMeshType, typename QuadratureType, typename Term>
void FiniteElementMethodTimeStepping<BasisOnMeshType, QuadratureType, Term>::
initialize(double timeStepWidth)
{
  LOG(DEBUG) << "FiniteElementMethodTimeStepping::initialize(timeStepWidth=" << timeStepWidth << ")";
  
  PetscErrorCode ierr;
  PetscScalar scale=-1.0/timeStepWidth;
  
  this->setMassMatrix();
  ierr=MatScale(this->data_.massMatrix(),scale);
  
  //this->preComputeSystemMatrix1();
  this->preComputeSystemMatrix(this->data_.systemMatrix());

  if (this->outputWriterManager_.hasOutputWriters())
  {
    LOG(WARNING) << "You have specified output writers for a FiniteElementMethod which is used for a time stepping problem. "
      "The output will not contain any solution data, only geometry. Probably you want to get output from the time stepping scheme, then define the output writers there.";
  }
  
  scale=-timeStepWidth;
  ierr=MatScale(this->data_.massMatrix(),scale); //required for explicit Euler but should be out commented for implicit1
}

/*
template<typename BasisOnMeshType, typename QuadratureType, typename Term>
bool FiniteElementMethodTimeStepping<BasisOnMeshType, QuadratureType, Term>::
setInitialValues(Vec& initialValues)
{
  LOG(TRACE)<<"FiniteElementMethodTimeStepping::setInitialValues()";
  
  dof_no_t nUnknowns = this->data_.nUnknowns();
  Vec &rightHandSide = this->data_.rightHandSide().values();
  
  std::vector<double> values;
  PythonUtility::getOptionVector(this->specificSettings_, "initialValues", nUnknowns, values);
  
  #ifndef NDEBUG
  LOG(DEBUG) << "Read in rhs values from config:";
  std::stringstream s;
  for(auto value : values)
  {
    s << " " << value;
  }
  LOG(DEBUG) << s.str();
  #endif
  
  PetscUtility::setVector(values, rightHandSide);
  
  return false;
  
}
*/
  
template<typename BasisOnMeshType, typename QuadratureType, typename Term>
constexpr int FiniteElementMethodTimeStepping<BasisOnMeshType, QuadratureType, Term>::
nComponents()
{
  return 1;   // this may be different for structural mechanics
}

template<typename BasisOnMeshType, typename QuadratureType, typename Term>
void FiniteElementMethodTimeStepping<BasisOnMeshType, QuadratureType, Term>::
checkDimensions(Mat &matrix, Vec &input)
{
#ifndef NDEBUG
  int nRows, nColumns;
  MatGetSize(matrix, &nRows, &nColumns);
  int nEntries;
  VecGetSize(input, &nEntries);
  if (nColumns != nEntries)
  {
    LOG(ERROR) << "matrix dimension " << nRows << "x" << nColumns << " does not match input vector (" << nEntries << ")!";
  }
  assert(nColumns == nEntries);
#endif
}

template<typename BasisOnMeshType, typename QuadratureType, typename Term>
bool FiniteElementMethodTimeStepping<BasisOnMeshType, QuadratureType, Term>::
knowsMeshType()
{
  return true;
}

template<typename BasisOnMeshType, typename QuadratureType, typename Term>
std::shared_ptr<Mesh::Mesh> FiniteElementMethodTimeStepping<BasisOnMeshType, QuadratureType, Term>::
mesh()
{
  return FiniteElementMethodBase<BasisOnMeshType, QuadratureType, Term>::mesh();
}

} // namespace SpatialDiscretization