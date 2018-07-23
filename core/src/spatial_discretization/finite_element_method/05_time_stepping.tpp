#include "spatial_discretization/finite_element_method/05_time_stepping.h"

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
  this->data_.initialize();
  this->setStiffnessMatrix();
  this->setMassMatrix();
  this->data_.finalAssembly();

  if (this->outputWriterManager_.hasOutputWriters())
  {
    LOG(WARNING) << "You have specified output writers for a FiniteElementMethod which is used for a time stepping problem. "
      "The output will not contain any solution data, only geometry. Probably you want to get output from the time stepping scheme, then define the output writers there.";
  }
}

template<typename BasisOnMeshType, typename QuadratureType, typename Term>
constexpr int FiniteElementMethodTimeStepping<BasisOnMeshType, QuadratureType, Term>::
nComponents()
{
  return 1;   // this may be different for structural mechanics
}


template<typename BasisOnMeshType, typename QuadratureType, typename Term>
void FiniteElementMethodTimeStepping<BasisOnMeshType, QuadratureType, Term>::
checkDimensions(Mat &stiffnessMatrix, Vec &input)
{
#ifndef NDEBUG
  int nRows, nColumns;
  MatGetSize(stiffnessMatrix, &nRows, &nColumns);
  int nEntries;
  VecGetSize(input, &nEntries);
  if (nColumns != nEntries)
  {
    LOG(ERROR) << "Stiffness matrix dimension " << nRows << "x" << nColumns << " does not match input vector (" << nEntries << ")!";
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