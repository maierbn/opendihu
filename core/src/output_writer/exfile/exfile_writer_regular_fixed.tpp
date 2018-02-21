#include "output_writer/exfile/exfile_writer.h"


namespace OutputWriter
{   

//! write exnode file to given stream
template<int D, typename BasisFunctionType>
void ExfileWriter<BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>>::
outputExelem(std::ostream &stream, std::vector<std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType>>> fieldVariables)
{
  
}

//! write exnode file to given stream
template<int D, typename BasisFunctionType>
void ExfileWriter<BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>>::
outputExnode(std::ostream &stream, std::vector<std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType>>> fieldVariables)
{
  
}

  
};  //namespace