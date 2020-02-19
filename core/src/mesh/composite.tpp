#include "mesh/composite.h"

namespace Mesh
{

template<int D>
CompositeOfDimension<D>::
CompositeOfDimension(PythonConfig specificSettings) :
{
  // created sub meshes
  for (int i = 0; i < nSubmeshes; i++)
  {
    meshes_[i] = std::make_shared<StructuredDeformableOfDimension<D>>(PythonConfig(specificSettings,i));
  }
}

//! get the total number of elements in the local domain
template<int D>
element_no_t CompositeOfDimension<D>::
nElementsLocal() const
{
  return nElementsLocal_;
}

//! get the total number of elements in the global domain
template<int D>
global_no_t CompositeOfDimension<D>::
nElementsGlobal() const
{
  return nElementsGlobal_;
}

//! the number of non-ghost nodes stored in the local partition
template<int D>
node_no_t CompositeOfDimension<D>::
nNodesLocalWithoutGhosts() const
{
  node_no_t nNodesLocalWithoutGhosts = 0;

  // iterate over submeshes and count number
  for(std::shared_ptr<StructuredDeformableOfDimension<D>> &submesh : meshes_)
  {
    nNodesLocalWithoutGhosts += submesh->nNodesLocalWithoutGhosts();
  }
  return nNodesLocalWithoutGhosts;
}

//! number of nodes in the mesh stored in the local partition, this also includes ghost nodes
template<int D>
node_no_t CompositeOfDimension<D>::
nNodesLocalWithGhosts() const
{
  node_no_t nNodesLocalWithGhosts = 0;

  // iterate over submeshes and count number
  for(std::shared_ptr<StructuredDeformableOfDimension<D>> &submesh : meshes_)
  {
    nNodesLocalWithGhosts += submesh->nNodesLocalWithGhosts();
  }
  return nNodesLocalWithGhosts;
}
