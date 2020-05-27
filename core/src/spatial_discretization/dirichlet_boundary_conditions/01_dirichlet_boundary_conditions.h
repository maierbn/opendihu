#pragma once

#include <Python.h>  // has to be the first included header
#include <vector>
#include <memory>

#include "spatial_discretization/dirichlet_boundary_conditions/00_dirichlet_boundary_conditions_base.h"
#include "partition/partitioned_petsc_mat/partitioned_petsc_mat.h"

namespace SpatialDiscretization
{

/** A class that handles Dirichlet type boundary conditions.
 *  In the python settings, the Dirichlet BC's need to be specified also for ghost dofs, when inputMeshIsGlobal=False.
  */
template<typename FunctionSpaceType, int nComponents>
class DirichletDirichletBoundaryConditionsBase :
  public DirichletBoundaryConditionsBase<FunctionSpaceType, nComponents>
{
public:

  //! use constructor of base class
  using DirichletBoundaryConditionsBase<FunctionSpaceType, nComponents>::DirichletBoundaryConditionsBase;

  //! set the Dirichlet boundary condition dofs to the values
  void applyInVector(std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nComponents>> fieldVariable);
};

template<typename FunctionSpaceType, int nComponents>
class DirichletBoundaryConditions :
  public DirichletDirichletBoundaryConditionsBase<FunctionSpaceType,nComponents>
{
public:

  //! use constructor of base class
  using DirichletDirichletBoundaryConditionsBase<FunctionSpaceType,nComponents>::DirichletDirichletBoundaryConditionsBase;

  typedef std::array<double,nComponents> ValueType;

  //! Set the boundary conditions in system matrix systemMatrixWrite, i.e. zero rows and columns of Dirichlet BC dofs and set diagonal to 1.
  //! Store the cleared matrix values (from systemMatrixRead) in boundaryConditionsRightHandSideSummand such that they can be used for adjusting the rhs vector afterwards.
  //! @param systemMatrixRead This matrix is used to read matrix values. It can be the same as systemMatrixWrite.
  //! @param systemMatrixWrite This matrix is changed, i.e. zeroed rows and columns.
  //! @param systemMatrixAlreadySet: if this is true, then the systemMatrixWrite is not changed. This is useful in a timestepping scheme where the dirichlet BC dofs do not change over time.
  void applyInSystemMatrix(const std::shared_ptr<PartitionedPetscMat<FunctionSpaceType>> systemMatrixRead,
                           std::shared_ptr<PartitionedPetscMat<FunctionSpaceType>> systemMatrixWrite,
                           std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nComponents>> boundaryConditionsRightHandSideSummand,
                           bool systemMatrixAlreadySet=false);

  //! set the boundary conditions to the right hand side using the information in boundaryConditionsRightHandSideSummand that was obtained by applyInSystemMatrix
  //! if the argument boundaryConditionsRightHandSideSummand is also set as rightHandSide, only write prescribed values into rightHandSide, do not add anything
  void applyInRightHandSide(std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nComponents>> rightHandSide,
                            std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nComponents>> boundaryConditionsRightHandSideSummand);

  //! This methods changes the prescribed boundary condition values to the values of a field variable. The dofs where BC values are prescribed remain the same.
  //! this is needed for the option "updatePrescribedValuesFromSolution"
  void updatePrescribedValuesFromSolution(std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nComponents>> solution);

protected:

  //! fill auxiliary ghost element data structures
  void initializeGhostElements();

  //! update boundary condition values in ownGhostElements_ from foreignGhostElements_
  void updateOwnGhostElements();

  struct GhostElement
  {
    std::vector<global_no_t> nonBoundaryConditionDofsOfRankGlobalPetsc;    //< the non-BC dofs of this element, as global petsc no. that are owned by the rank with no neighbouringRankNo
    std::vector<global_no_t> boundaryConditionDofsGlobalPetsc;     //< the Dirichlet BC dofs of this element
    std::vector<ValueType> boundaryConditionValues;                //< the prescribed value, corresponding to boundaryConditionDofsGlobalPetsc
  };
  
  std::map<int,std::vector<GhostElement>> foreignGhostElements_;   //< ghost elements that are normal elements on this rank, key is the rankNo of the rank to send them to
  std::vector<GhostElement> ownGhostElements_;                     //< the ghost elements for this rank

  std::vector<std::pair<int,int>> nElementsFromRanks_;             //< helper variable, (foreignRank,nElements), number of elements to receive from foreignRank, will be initialized by initializeGhostElements() and used by updateOwnGhostElements().
};

} // namespace

#include "spatial_discretization/dirichlet_boundary_conditions/01_dirichlet_boundary_conditions.tpp"
