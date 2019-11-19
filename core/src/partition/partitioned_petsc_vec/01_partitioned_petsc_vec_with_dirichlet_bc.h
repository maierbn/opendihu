#pragma once

#include <memory>

#include "control/types.h"
#include "partition/rank_subset.h"
#include "function_space/function_space.h"
#include "mesh/type_traits.h"
#include "partition/partitioned_petsc_vec/partitioned_petsc_vec_base.h"
#include "spatial_discretization/boundary_conditions/dirichlet_boundary_conditions.h"

// forward declaration
namespace FunctionSpace
{
template<typename MeshType, typename BasisFunctionType>
class FunctionSpace;
}

/** This vector accounts for prescribed Dirichlet BC values. Only non-BC values are stored in the vector.
 *  The getValues and setValues methods handle normal indices, without special treatment of BC dofs.
 *  Internally, setting values for Dirichlet BC dofs has no effect and getting values for Dirichlet BC dofs gives the actual
 *  prescribed values. The internal data representation, (representationCombinedLocal and representationCombinedGlobal) uses
 *  a new ordering that simply skips Dirichlet BC dofs. Furthermore, the values for different components are stored after each other:
 *  Example: 5 dofs, 2 components, dof2c0 and dof4c1 have prescribed Dirichlet values, then the vector will contain
 *    [dof0c0, dof1c0, dof3c0, dof4c0, dof0c1, dof1c1, dof2c1, dof3c1] in the serial case, plus ghost values in the parallel case.
 * The nComponentsDirichletBc template parameter has no special meaning, it is only used for the derived class PartitionedPetscVecForHyperelasticity
 */
template<typename FunctionSpaceType, int nComponents, int nComponentsDirichletBc=nComponents>
class PartitionedPetscVecWithDirichletBc:
  public PartitionedPetscVec<FunctionSpaceType,nComponents>
{
public:
 
  //! constructor
  PartitionedPetscVecWithDirichletBc(std::shared_ptr<Partition::MeshPartition<FunctionSpaceType>> meshPartition,
                                     std::shared_ptr<SpatialDiscretization::DirichletBoundaryConditions<FunctionSpaceType,nComponentsDirichletBc>> dirichletBoundaryConditions,
                                     std::string name);
 
  //! Communicates the ghost values from the global vectors to the local vector and sets the representation to local.
  //! The representation has to be global, afterwards it is set to local.
  void startGhostManipulation();
  
  //! Communicates the ghost values from the local vectors back to the global vector and sets the representation to global.
  //! The representation has to be local, afterwards it is set to global.
  void finishGhostManipulation();
  
  //! zero all values in the local ghost buffer. Needed if between startGhostManipulation() and finishGhostManipulation() only some ghost will be reassigned. To prevent that the "old" ghost values that were present in the local ghost values buffer get again added to the real values which actually did not change.
  void zeroGhostBuffer();

  //! wrapper to the PETSc VecSetValues, acting on the local data or global data, the indices ix are the local dof nos
  void setValues(int componentNo, PetscInt ni, const PetscInt ix[], const PetscScalar y[], InsertMode iora = INSERT_VALUES);
  
  //! wrapper to the PETSc VecSetValue, acting on the local data or global data, the row is local dof no
  void setValue(int componentNo, PetscInt row, PetscScalar value, InsertMode mode = INSERT_VALUES);

  //! wrapper to the PETSc VecGetValues, acting on the local data or global data, the indices ix are the local dof nos
  void getValues(int componentNo, PetscInt ni, const PetscInt ix[], PetscScalar y[]) const;

  //! get values from a different vector but using the indexing for the global vector
  void getValuesGlobal(Vec vector, int componentNo, PetscInt ni, const PetscInt ix[], PetscScalar y[]) const;

  //! get a single value
  double getValue(int componentNo, PetscInt row) const;

  //! set all entries to zero, wraps VecZeroEntries
  void zeroEntries();
  
  //! output the vector to stream, for debugging
  void output(std::ostream &stream) const;

  //! get the index in the internal vector in global numbering, from componentNo and local dof no
  global_no_t nonBCDofNoGlobal(int componentNo, dof_no_t localDofNo) const;

  //! get number of local entries
  int nEntriesLocal() const;

  //! get number of global entries
  global_no_t nEntriesGlobal() const;

  //! get the Petsc Vec that contains all components but no values for dirichlet BC dofs
  Vec &valuesGlobal();

  //! get the index in the internal vector in local numbering, from componentNo and local dof no
  dof_no_t nonBCDofNoLocal(int componentNo, dof_no_t localDofNo) const;

  //! determine if the dof entry is prescribed
  bool isPrescribed(int componentNo, dof_no_t localDofNo) const;

  //! get a reference to the internal dofNoLocalToDofNoNonBcGlobal_ data structure
  const std::array<std::vector<dof_no_t>,nComponents> &dofNoLocalToDofNoNonBcGlobal() const;

  //! set the internal representation to be combined global, i.e. using the global vector (vectorCombinedWithoutDirichletDofsGlobal_), if it was local, ghost buffer entries are discarded (use finishGhostManipulation to consider ghost dofs)
  void setRepresentationGlobal();

protected:

  //! prepare internal variables such that the vector can be created by createVector afterwards
  //! @param nComponents_ the actual number of components to initialize, this is usually the number of components of the whole vector, it is just different for PartitionedPetscVecForHyperelasticity
  //! @param offsetInGlobalNumberingPerRank Unused numbers in the non-bc global numbering at the end of the local numbers of each rank, only needed for PartitionedPetscVecForHyperelasticity
  void initialize(int nComponents_, int offsetInGlobalNumberingPerRank = 0);

  //! create the distributed Petsc vector vectorCombinedWithoutDirichletDofs_
  void createVector();

  //! set the internal representation to be combined local, i.e. using the local vector (vectorCombinedWithoutDirichletDofsLocal_), ghost buffer is not filled (use startGhostManipulation to consider ghost dofs)
  void setRepresentationLocal();

  std::shared_ptr<SpatialDiscretization::DirichletBoundaryConditions<FunctionSpaceType,nComponentsDirichletBc>> dirichletBoundaryConditions_; //< the dirichlet boundary conditions object that contains dofs and values of Dirichlet BCs

  /** The local vector contains the nodal/dof values for the local portion of the current rank. This includes ghost nodes.
   *  The global vector manages the whole data and also only stores the local portion of it on the current rank.
   *  It shares the memory with the local vector. The local vector uses local indices whereas the global vector is accessed using globalPetsc indexing.
   *  If we want to manipulate data, we have to ensure consistency in the parallel global vector (VecGhostUpdateBegin,VecGhostUpdateEnd)
   *  and then fetch the local portion of the global vector together with ghost values into a local vector (VecGhostGetLocalForm),
   *  then manipulate the values (using standard VecSetValues/VecGetValues routines) and commit the results (VecGhostRestoreLocalForm, VecGhostUpdateBegin, VecGhostUpdateEnd).
   */
  Vec vectorCombinedWithoutDirichletDofsGlobal_;  ///< the values of all components in "struct of array" ordering, Dirichlet BC dofs are left out. This is accessed using global ordering.
  Vec vectorCombinedWithoutDirichletDofsLocal_;  ///< the values of all components in "struct of array" ordering, Dirichlet BC dofs are left out. This is accessed using local ordering.

  std::array<int,nComponents> nNonBcDofsWithoutGhosts_;  ///< the local without ghosts number of entries in the vector, without the Dirichlet BC dofs
  int nNonBcDofsGhosts_;         ///< number of ghost values
  int nEntriesLocal_;    ///< the local number of entries in the vector, without Ghost dofs
  int nEntriesGlobal_;  ///< the global number of entries in the vector, without the Dirichlet BC dofs
  global_no_t nonBcDofNoGlobalBegin_;   ///< the first no in the non-bc global numbering
  int nDofsLocal_;     ///< same as nEntriesLocal_ (needed for PartitionedPetscVecForHyperelasticity)

  std::vector<int> nonBcGhostDofNosGlobal_;  ///< non-bc ghost dofs in non-bc global indexing

  std::array<std::vector<dof_no_t>,nComponents> dofNoLocalToDofNoNonBcGlobal_;   ///< mapping from component no and local dof no to the numbering used for the combined vector, for local dofs with ghosts
  std::array<std::vector<dof_no_t>,nComponents> dofNoLocalToDofNoNonBcLocal_;    ///< mapping from component no and local dof no to the local number of the non-bc dof numbering
  std::array<std::vector<double>,nComponents> boundaryConditionValues_;   ///< prescribed boundary condition values for local dof nos (normal local dof numbering)
  std::array<std::vector<bool>,nComponents> isPrescribed_;                ///< for every local dof no, if the dof has a prescribed Dirichlet BC value
};

template<typename FunctionSpaceType, int nComponents>
std::ostream &operator<<(std::ostream &stream, const PartitionedPetscVecWithDirichletBc<FunctionSpaceType,nComponents> &vector);

#include "partition/partitioned_petsc_vec/01_partitioned_petsc_vec_with_dirichlet_bc.tpp"
