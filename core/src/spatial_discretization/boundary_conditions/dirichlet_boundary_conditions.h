#pragma once

#include <Python.h>  // has to be the first included header
#include <vector>
#include <memory>

#include "spatial_discretization/boundary_conditions/boundary_conditions_base.h"
#include "partition/partitioned_petsc_mat/partitioned_petsc_mat.h"

namespace SpatialDiscretization
{

/** A class that handles Dirichlet type boundary conditions.
 *  In the python settings, the Dirichlet BC's need to be specified also for ghost dofs, when inputMeshIsGlobal=False.
  */
template<typename FunctionSpaceType, int nComponents>
class DirichletBoundaryConditionsBase :
  public BoundaryConditionsBase<FunctionSpaceType, nComponents>
{
public:

  //! use constructor of base class
  using BoundaryConditionsBase<FunctionSpaceType, nComponents>::BoundaryConditionsBase;

  //! set the Dirichlet boundary condition dofs to the values
  void applyInVector(std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nComponents>> fieldVariable);
};

template<typename FunctionSpaceType, int nComponents>
class DirichletBoundaryConditions :
  public DirichletBoundaryConditionsBase<FunctionSpaceType,nComponents>
{
public:

  //! use constructor of base class
  using DirichletBoundaryConditionsBase<FunctionSpaceType,nComponents>::DirichletBoundaryConditionsBase;

  typedef std::array<double,nComponents> ValueType;

  //! set the boundary conditions to system matrix, i.e. zero rows and columns of Dirichlet BC dofs and set diagonal to 1. Store the cleared matrix values in boundaryConditionsRightHandSideSummand such that they can be used for adjusting the rhs vector afterwards
  //! @param systemMatrixAlreadySet: if this is true, then the systemMatrix is not changed. This is useful in a timestepping scheme where the dirichlet BC dofs do not change over time.
  void applyInSystemMatrix(std::shared_ptr<PartitionedPetscMat<FunctionSpaceType>> systemMatrix,
                           std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nComponents>> boundaryConditionsRightHandSideSummand,
                           bool systemMatrixAlreadySet=false);

  //! set the boundary conditions to the right hand side using the information in boundaryConditionsRightHandSideSummand that was obtained by applyInSystemMatrix
  //! if the argument boundaryConditionsRightHandSideSummand is also set as rightHandSide, only write prescribed values into rightHandSide, do not add anything
  void applyInRightHandSide(std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nComponents>> rightHandSide,
                            std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nComponents>> boundaryConditionsRightHandSideSummand);

protected:

  //! fill auxiliary ghost element data structures
  void initializeGhostElements();

  struct GhostElement
  {
    std::vector<global_no_t> nonBoundaryConditionDofsOfRankGlobalPetsc;    ///< the non-BC dofs of this element, as global petsc no. that are owned by the rank with no neighbouringRankNo
    std::vector<global_no_t> boundaryConditionDofsGlobalPetsc;      ///< the Dirichlet BC dofs of this element
    std::vector<ValueType> boundaryConditionValues;   ///< the prescribed value, corresponding to boundaryConditionDofsGlobalPetsc
  };
  
  std::map<int,std::vector<GhostElement>> foreignGhostElements_;   ///< ghost elements that are normal elements on this rank, key is the rankNo of the rank to send them to
  std::vector<GhostElement> ownGhostElements_;   ///< the ghost elements for this rank
};

} // namespace

#include "spatial_discretization/boundary_conditions/dirichlet_boundary_conditions.tpp"
