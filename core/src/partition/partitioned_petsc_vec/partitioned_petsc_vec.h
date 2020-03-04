#pragma once

#include <memory>

#include "control/types.h"
#include "partition/rank_subset.h"
#include "function_space/function_space.h"
#include "mesh/type_traits.h"
#include "partition/partitioned_petsc_vec/partitioned_petsc_vec_base.h"

// forward declaration
namespace FunctionSpace
{
template<typename MeshType, typename BasisFunctionType>
class FunctionSpace;
}

/** This encapsulates a Petsc Vec, combined with the partition of the mesh.
 *  For each component a local Vec is stored that holds all values of that component.
 *  Global Petsc numbering: such that each rank has its own contiguous subset in this numbering of the total range.
 *  Global Natural numbering: normal indexing proceeding fastest in x, then in y, then in z direction, over the whole domain.
 *  Local numbering: starting with 0, first all non-ghost values, then the ghost indices.
 * *
 *  This particular standard specialization is for non-structured meshes or no meshes and currently completely serial, 
 *  it is the placeholder as long as the partial specialization for unstructured meshes is not implemented. (It will never be)
 *  This means some of the methods here have no effect. Most likely you want to have a look at the partial specialization for
 *  structured meshes, below this class.
 */
template<typename FunctionSpaceType, int nComponents, typename = typename FunctionSpaceType::Mesh>
class PartitionedPetscVec : 
  public PartitionedPetscVecBase<FunctionSpaceType>
{
public:
  //! constructor
  PartitionedPetscVec(std::shared_ptr<Partition::MeshPartition<FunctionSpaceType,typename FunctionSpaceType::Mesh>> meshPartition, std::string name);
 
  //! constructor, copy values from existing rhs vector or reuse Petsc Vec's from rhs vector
  //! @param reuseData if true, it uses the same Petsc Vec's for global and local data. If false, it creates new Petsc Vec's and copies the values.
  template<int nComponents2>
  PartitionedPetscVec(PartitionedPetscVec<FunctionSpaceType,nComponents2> &rhs, std::string name, bool reuseData=false);
  
  //! this has to be called before the vector is manipulated (i.e. VecSetValues or vecZeroEntries is called)
  void startGhostManipulation();
  
  //! this has to be called after the vector is manipulated (i.e. VecSetValues or vecZeroEntries is called)
  void finishGhostManipulation();
  
  //! zero all values in the local ghost buffer. Needed if between startGhostManipulation() and finishGhostManipulation() only some ghost will be reassigned. To prevent that the "old" ghost values that were present in the local ghost values buffer get again added to the real values which actually did not change.
  void zeroGhostBuffer();

  //! wrapper to the PETSc VecSetValues, acting only on the local data, the indices ix are the local dof nos
  void setValues(int componentNo, PetscInt ni, const PetscInt ix[], const PetscScalar y[], InsertMode iora);
  
  //! wrapper to the PETSc VecSetValue, acting only on the local data
  void setValue(int componentNo, PetscInt row, PetscScalar value, InsertMode mode);
  
  //! set values from another vector, only the first components are copied, if nComponents != nComponents2
  template<int nComponents2>
  void setValues(PartitionedPetscVec<FunctionSpaceType,nComponents2> &rhs);

  //! wrapper to the PETSc VecGetValues, acting only on the local data, the indices ix are the local dof nos
  void getValues(int componentNo, PetscInt ni, const PetscInt ix[], PetscScalar y[]);
  
  //! wrapper to the PETSc VecGetValues, on the global vector with global/Petsc indexing. This is not the global natural numbering!
  void getValuesGlobalPetscIndexing(int componentNo, PetscInt ni, const PetscInt ix[], PetscScalar y[]);
  
  //! set all entries to zero, wraps VecZeroEntries
  void zeroEntries();
  
  //! get the internal PETSc vector values, the local vector for the specified component
  Vec &valuesLocal(int componentNo = 0);

  //! get the internal PETSc vector values, the global vector for the specified component
  Vec &valuesGlobal(int componentNo);

  //! if the vector has multiple components, return a nested Vec of the global vector, else return the global vector
  Vec &valuesGlobal();

  //! fill a contiguous vector with all components after each other, "struct of array"-type data layout.
  //! after manipulation of the vector has finished one has to call restoreValuesContiguous
  Vec &getValuesContiguous();

  //! copy the values back from a contiguous representation where all components are in one vector to the standard internal format of PartitionedPetscVec where there is one local vector with ghosts for each component.
  //! this has to be called
  void restoreValuesContiguous();

  //! set the internal representation to be global, for unstructured meshes this means "not contiguous", because there is no local or global vector
  void setRepresentationGlobal();

  //! set the internal representation to be local, for unstructured meshes this means "not contiguous", because there is no local or global vector
  void setRepresentationLocal();

  //! set the internal representation to be contiguous, i.e. using the contiguous vectors
  void setRepresentationContiguous();

  //! output the vector to stream, for debugging
  void output(std::ostream &stream);

  template<typename FunctionSpaceType2, int nComponents2, typename Dummy> friend class PartitionedPetscVec;

protected:
  
  //! create the values vectors
  void createVector();
  
  std::array<Vec,nComponents> values_;  ///< the (serial) Petsc vectors that contains all the data, one for each component
  Vec valuesContiguous_ = PETSC_NULL;   ///< global vector that has all values of the components concatenated, i.e. in a "struct of arrays" memory layout
  Vec vectorNestedGlobal_;       ///< a VecNest object containing the global values, only in used if nComponents > 1
};

/** This is the partial specialization of Petsc Vec's on a mesh for structured meshes.
 *  This class represents the values of a field which is defined on the domain discretized by a mesh and basis function type.
 *  The field can be scalar or have multiple components.
 *  This means that for every degree of freedom (=node for Lagrange ansatz functions) its stores as many values as components.
 *
 *  This class wraps the normal Petsc Vec and mainly wraps the setValues and getValues functions.
 *  There is one Vec per component.
 *  There is a local and a global version of the vecs, vectorLocal_ and vectorGlobal_. The local vectors only store the values
 *  that are on the local domain. The global vector also only stores the local values but additionally the ghost values.
 *  The global vector is the one that should be used in any computation with Petsc functions on the data. The local vectors
 *  will internally be used when setValues and getValues are called with local indices, whereas the global vector needs global indices.
 *  Petsc does the data transfer between the local and the global vector and fills and accumulates the correct ghost values.
 *
 *  Only one of the local or global vector contain the currently valid data at any time, because they even share memory.
 *  Which one this is is save in this->currentRepresentation_ (read the comments there). The transformations between representations
 *  is done internally and should not be cared about from the outside of this class
 *  (except you have to properly call startGhostManipulation() and finishGhostManipulation()).
 *
 *  valuesGlobal(componentNo) returns the global Petsc Vec of the given component. valuesGlobal() returns a nested Petsc Vec of all components.
 *
 *  Internally, an own DMDA object is generated, separately from the one in MeshPartition. This object now refers to nodes (as opposite the one of MeshPartition which refers to elements).
 *  This object is created such that it matches the partition given by the meshPartition.
 */
template<typename MeshType, typename BasisFunctionType, int nComponents>
class PartitionedPetscVecNComponentsStructured : 
  public PartitionedPetscVecBase<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>>
{
public:
 
  //! constructor, construct a petsc Vec with meshPartition that can hold the values for a field variable with nComponents components
  PartitionedPetscVecNComponentsStructured(std::shared_ptr<Partition::MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,MeshType>> meshPartition, std::string name);
 
  //! constructor, copy values from existing rhs vector or reuse Petsc Vec's from rhs vector
  //! @param reuseData if true, it uses the same Petsc Vec's for global and local data. If false, it creates new Petsc Vec's and copies the values.
  //! reuseData should be used with care, because the internal representation will also be copied
  template<int nComponents2>
  PartitionedPetscVecNComponentsStructured(PartitionedPetscVec<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,nComponents2> &rhs, std::string name, bool reuseData=false);
 
  //! Communicates the ghost values from the global vectors to the local vector and sets the representation to local.
  //! The representation has to be global, afterwards it is set to local.
  void startGhostManipulation();
  
  //! Communicates the ghost values from the local vectors back to the global vector and sets the representation to global.
  //! The representation has to be local, afterwards it is set to global.
  void finishGhostManipulation();
  
  //! zero all values in the local ghost buffer. Needed if between startGhostManipulation() and finishGhostManipulation() only some ghost will be reassigned. To prevent that the "old" ghost values that were present in the local ghost values buffer get again added to the real values which actually did not change.
  void zeroGhostBuffer();

  //! wrapper to the PETSc VecSetValues, acting only on the local data, the indices ix are the local dof nos
  void setValues(int componentNo, PetscInt ni, const PetscInt ix[], const PetscScalar y[], InsertMode iora);
  
  //! wrapper to the PETSc VecSetValue, acting only on the local data
  void setValue(int componentNo, PetscInt row, PetscScalar value, InsertMode mode);

  //! for a single component vector set all values. They have to be enough for all local dof including ghosts.
  void setValuesWithGhosts(int componentNo, std::vector<double> &values, InsertMode petscInsertMode);
  
  //! for a single component vector set all values. values does not contain ghost dofs.
  void setValuesWithoutGhosts(int componentNo, std::vector<double> &values, InsertMode petscInsertMode);

  //! set values from another vector, only the first components are copied, if nComponents != nComponents2
  template<int nComponents2>
  void setValues(PartitionedPetscVec<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,nComponents2> &rhs);

  //! set values of a specific component from another vector, this is the opposite operation to extractComponentCopy
  void setValues(int componentNo, std::shared_ptr<PartitionedPetscVec<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,1>> fieldVariable);

  //! set the values for the given component from a petsc Vec, name is only for debugging output
  void setValues(int componentNo, Vec petscVector, std::string name = "");

  //! extract a single component, this field variable can have any representation
  //! It set the representation of extractedPartitionedPetscVec to local.
  void extractComponentCopy(int componentNo, std::shared_ptr<PartitionedPetscVec<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,1>> extractedPartitionedPetscVec);

  //! extract a component from the shared vector (no copy), this field variable cannot be used any longer and is set to invalid, until restoreExtractedComponent is called.
  void extractComponentShared(int componentNo, std::shared_ptr<PartitionedPetscVec<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,1>> extractedPartitionedPetscVec);

  //! restore the extracted raw array to petsc and make the field variable usable again
  template<int nComponents2>
  void restoreExtractedComponent(std::shared_ptr<PartitionedPetscVec<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,nComponents2>> extractedPartitionedPetscVec);

  //! wrapper to the PETSc VecGetValues, acting only on the local data, the indices ix are the local dof nos
  void getValues(int componentNo, PetscInt ni, const PetscInt ix[], PetscScalar y[]);
  
  //! wrapper to the PETSc VecGetValues, on the global vector with global/Petsc indexing. This is not the global natural numbering!
  void getValuesGlobalPetscIndexing(int componentNo, PetscInt ni, const PetscInt ix[], PetscScalar y[]);
  
  //! get all locally stored values, i.e. with ghosts
  void getLocalValues(int componentNo, std::vector<double> &values);
  
  //! set all entries to zero, wraps VecZeroEntries
  void zeroEntries();
  
  //! get the local Vector of a specified component
  Vec &valuesLocal(int componentNo = 0);
  
  //! get the global Vector of a specified component
  Vec &valuesGlobal(int componentNo);
  
  //! if the vector has multiple components, return a nested Vec of the global vector, else return the global vector
  Vec &valuesGlobal();

  //! fill a contiguous vector with all components after each other, "struct of array"-type data layout.
  //! after manipulation of the vector has finished one has to call restoreValuesContiguous
  virtual Vec &getValuesContiguous() = 0;

  //! copy the values back from a contiguous representation where all components are in one vector to the standard internal format of PartitionedPetscVec where there is one local vector with ghosts for each component.
  //! this has to be called
  virtual void restoreValuesContiguous() = 0;

  //! set the internal representation to be global, i.e. using the global vectors, if it was local, ghost buffer entries are discarded (use finishGhostManipulation to consider ghost dofs)
  void setRepresentationGlobal();

  //! set the internal representation to be local, i.e. using the local vectors, ghost buffer is not filled (use startGhostManipulation to consider ghost dofs)
  void setRepresentationLocal();

  //! set the internal representation to be contiguous, i.e. using the contiguous vectors
  void setRepresentationContiguous();

  //! get a vector of local dof nos (from meshPartition), without ghost dofs
  std::vector<PetscInt> &localDofNosWithoutGhosts();
  
  //! output the vector to stream, for debugging
  void output(std::ostream &stream);

  template<typename FunctionSpaceType2, int nComponents2, typename Dummy> friend class PartitionedPetscVec;
  template<typename MeshType2, typename BasisFunctionType2, int nComponents2> friend class PartitionedPetscVecNComponentsStructured;

protected:
 
  //! create a distributed Petsc vector, according to partition
  void createVector();
  
  std::shared_ptr<DM> dm_;    ///< PETSc DMDA object that stores topology information and everything needed for communication of ghost values
  bool ghostManipulationStarted_;   ///< if startGhostManipulation() was called but not yet finishGhostManipulation(). This indicates that finishGhostManipulation() can be called next without giving an error.
  
  std::array<Vec,nComponents> vectorLocal_;   ///< local vector that holds the local Vecs, is filled by startGhostManipulation and can the be manipulated, afterwards the results need to get copied back by finishGhostManipulation
  std::array<Vec,nComponents> vectorGlobal_;  ///< the global distributed vector that holds the actual data
  Vec valuesContiguous_ = PETSC_NULL;   ///< global vector that has all values of the components concatenated, i.e. in a "struct of arrays" memory layout. This is never used if nComponents = 1

  std::vector<PetscInt> temporaryIndicesVector_;   ///< a temporary vector that will be used whenever indices are to be computed, this avoids creating and deleting local vectors which is time-consuming (found out by perftools on hazelhen)

  const double *extractedData_ = nullptr;   ///< the data array of valuesContiguous_, used when a component is extracted by extractComponentShared, then the representation is set to invalid
  std::vector<double> savedValues_;   ///< temporary storage of values that would be overwritten by ghost value operations of the extracted field variable
  Vec savedVectorLocal_;        ///< when this PartitionedPetscVec has nComponents=1 and extractComponentShared is called, there is no valuesContiguous_ vector in use (because it is only one component anyway, replacement is globalVector_[0]). Then the extracted field variable gets copies of the own vectorLocal_ and vectorGlobal_ set, the original pointer vectorLocal_ and vectorGlobal_ are saved in this variable and reset when restoreValuesContiguous is called.
  Vec savedVectorGlobal_;        ///< when this PartitionedPetscVec has nComponents=1 and extractComponentShared is called, there is no valuesContiguous_ vector in use (because it is only one component anyway, replacement is globalVector_[0]). Then the extracted field variable gets copies of the own vectorLocal_ and vectorGlobal_ set, the original pointer vectorLocal_ and vectorGlobal_ are saved in this variable and reset when restoreValuesContiguous is called.
  Vec vectorNestedGlobal_;       ///< a VecNest object containing the global values, only in used if nComponents > 1
};

/** This is a partial specialization for structured meshes with multiple components.
 */
template<typename MeshType, typename BasisFunctionType, int nComponents>
class PartitionedPetscVec<
  FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,
  nComponents,
  Mesh::isStructuredOrComposite<MeshType>> :
  public PartitionedPetscVecNComponentsStructured<MeshType,BasisFunctionType,nComponents>
{
public:
  using PartitionedPetscVecNComponentsStructured<MeshType,BasisFunctionType,nComponents>::PartitionedPetscVecNComponentsStructured;
  
  //! fill a contiguous vector with all components after each other, "struct of array"-type data layout.
  //! after manipulation of the vector has finished one has to call restoreValuesContiguous
  Vec &getValuesContiguous() override;
  
  //! copy the values back from a contiguous representation where all components are in one vector to the standard internal format of PartitionedPetscVec where there is one local vector with ghosts for each component.
  //! this has to be called
  void restoreValuesContiguous() override;
};

/** This is the partial specialization for 1-component vectors
 */
template<typename MeshType, typename BasisFunctionType>
class PartitionedPetscVec<
  FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,
  1,
  Mesh::isStructuredOrComposite<MeshType>> :
  public PartitionedPetscVecNComponentsStructured<MeshType,BasisFunctionType,1>
{
public:
  using PartitionedPetscVecNComponentsStructured<MeshType,BasisFunctionType,1>::PartitionedPetscVecNComponentsStructured;
  
  //! fill a contiguous vector with all components after each other, "struct of array"-type data layout.
  //! after manipulation of the vector has finished one has to call restoreValuesContiguous
  Vec &getValuesContiguous() override;
  
  //! copy the values back from a contiguous representation where all components are in one vector to the standard internal format of PartitionedPetscVec where there is one local vector with ghosts for each component.
  //! this has to be called
  void restoreValuesContiguous() override;
};

template<typename FunctionSpaceType, int nComponents>
std::ostream &operator<<(std::ostream &stream, PartitionedPetscVec<FunctionSpaceType,nComponents> &vector);

#include "partition/partitioned_petsc_vec/partitioned_petsc_vec_unstructured.tpp"
#include "partition/partitioned_petsc_vec/partitioned_petsc_vec_n_components_structured.tpp"
#include "partition/partitioned_petsc_vec/partitioned_petsc_vec_structured.tpp"
