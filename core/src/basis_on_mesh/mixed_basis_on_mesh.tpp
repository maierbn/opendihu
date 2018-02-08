#include "basis_on_mesh/mixed_basis_on_mesh.h"

#include <Python.h>  // has to be the first included header

namespace BasisOnMesh
{

template<typename LowOrderBasisOnMeshType,typename HighOrderBasisOnMeshType>
Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>::
Mixed(PyObject *specificSettings) : Mesh::Mesh(specificSettings), 
  lowOrderBasisOnMesh_(std::make_shared<LowOrderBasisOnMeshType>(specificSettings, true)), 
  highOrderBasisOnMesh_(std::make_shared<HighOrderBasisOnMeshType>(specificSettings, true))
{
}
  
template<typename LowOrderBasisOnMeshType,typename HighOrderBasisOnMeshType>
void Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>::
initialize()
{
  lowOrderBasisOnMesh_->initialize();
  highOrderBasisOnMesh_->initialize();
}
  
template<typename LowOrderBasisOnMeshType,typename HighOrderBasisOnMeshType>
int Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>::
dimension() const
{
  return highOrderBasisOnMesh_->dimension();
}
  
template<typename LowOrderBasisOnMeshType,typename HighOrderBasisOnMeshType>
constexpr int Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>::
dim()
{
  return HighOrderBasisOnMeshType::Mesh::dim();
}
  
template<typename LowOrderBasisOnMeshType,typename HighOrderBasisOnMeshType>
node_no_t Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>::
nNodes() const
{
  return highOrderBasisOnMesh_->nNodes();
}

template<typename LowOrderBasisOnMeshType,typename HighOrderBasisOnMeshType>
std::shared_ptr<LowOrderBasisOnMeshType> Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>::
lowOrderBasisOnMesh()
{
  return lowOrderBasisOnMesh_;
}

template<typename LowOrderBasisOnMeshType,typename HighOrderBasisOnMeshType>
std::shared_ptr<HighOrderBasisOnMeshType> Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>::
highOrderBasisOnMesh()
{
  return highOrderBasisOnMesh_;
}

};  //namespace