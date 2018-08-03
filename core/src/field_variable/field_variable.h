#pragma once

#include "field_variable/08_field_variable_vector.h"

namespace FieldVariable
{

/** General field variable
 */
template<typename BasisOnMeshType,int nComponents>
class FieldVariable :
  public FieldVariableVector<BasisOnMeshType,nComponents>
{
public:
  //! inherited constructor
  using FieldVariableVector<BasisOnMeshType,nComponents>::FieldVariableVector;

  typedef BasisOnMeshType BasisOnMesh;
  
  //! this has to be called before the vector is manipulated (i.e. VecSetValues or vecZeroEntries is called), to ensure that the current state of the vector is fetched from the global vector
  void startVectorManipulation();
  
  //! this has to be called after the vector is manipulated (i.e. VecSetValues or vecZeroEntries is called), to ensure that operations on different partitions are merged by Petsc
  void finishVectorManipulation();
};


// output operator
template<typename BasisOnMeshType,int nComponents>
std::ostream &operator<<(std::ostream &stream, const FieldVariable<BasisOnMeshType,nComponents> &rhs)
{
  rhs.output(stream);
  return stream;
}

};  // namespace

#include "field_variable/field_variable.tpp"