#pragma once

#include "field_variable/09_field_variable_composite.h"

namespace FieldVariable
{

/** General field variable, this is a vector with as many entries as there are unknowns in the function space.
 *  It uses a PartitionedPetscVec at data container which is a wrapper for Petsc Vec.
 */
template<typename FunctionSpaceType,int nComponents>
class FieldVariable :
  public FieldVariableComposite<FunctionSpaceType,nComponents>
{
public:
  //! inherited constructor
  using FieldVariableComposite<FunctionSpaceType,nComponents>::FieldVariableComposite;

  typedef FunctionSpaceType FunctionSpace;
  
  //! this has to be called before the vector is manipulated (i.e. VecSetValues or vecZeroEntries is called), to ensure that the current state of the vector is fetched from the global vector
  void startGhostManipulation();
  
  //! zero all values in the local ghost buffer. Needed if between startGhostManipulation() and finishGhostManipulation() only some ghost will be reassigned. To prevent that the "old" ghost values that were present in the local ghost values buffer get again added to the real values which actually did not change.
  void zeroGhostBuffer();

  //! this has to be called after the vector is manipulated (i.e. VecSetValues or vecZeroEntries is called), to ensure that operations on different partitions are merged by Petsc
  //! It sums up the values in the ghost buffer and the actual nodal value.
  void finishGhostManipulation();

  //! set the internal representation to be global, i.e. using the global vectors, if it was local, ghost buffer entries are discarded (use finishGhostManipulation to consider ghost dofs)
  void setRepresentationGlobal();

  //! set the internal representation to be local, i.e. using the local vectors, ghost buffer is not filled (use startGhostManipulation to consider ghost dofs)
  void setRepresentationLocal();

  //! set the internal representation to be contiguous, i.e. using the contiguous vectors
  void setRepresentationContiguous();

  //! get the current internal data representation, or `noVector` if the internal vector is not initialized.
  Partition::values_representation_t currentRepresentation() const;

  //! set the internal representation. Allows to reset the representation without writing back values if they are unchanged since the vector was in state `representation` the last time.
  void setRepresentation(Partition::values_representation_t representation, values_modified_t values = values_modified_t::values_modified);

  //! check if the field variable contains Nan or Inf values
  bool containsNanOrInf();
};

// output operator
template<typename FunctionSpaceType,int nComponents>
std::ostream &operator<<(std::ostream &stream, const FieldVariable<FunctionSpaceType,nComponents> &rhs)
{
#ifndef NDEBUG
  rhs.output(stream);
#endif
  return stream;
}

} // namespace

#include "field_variable/field_variable.tpp"
