#include "output_writer/paraview/get_connectivity_values_unstructured_mesh.h"

namespace OutputWriter
{

template<typename FunctionSpaceType>
void GetConnectivityValuesUnstructuredMesh<FunctionSpaceType>::
get(std::shared_ptr<FunctionSpaceType> functionSpace, std::vector<int> &connectivityValues)
{
  // do not do anything for not UnstructuredDeformableOfDimension mesh
}

template<int D, typename BasisFunctionType>
void GetConnectivityValuesUnstructuredMesh<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>>::
get(std::shared_ptr<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>> functionSpace, std::vector<int> &connectivityValues)
{
  LOG(DEBUG) << "GetConnectivityValuesUnstructuredMesh, " << functionSpace->nElementsLocal() << " elements";
  for (element_no_t elementNo = 0; elementNo < functionSpace->nElementsLocal(); elementNo++)
  {
    //const std::vector<int> &nodeGlobalNo = functionSpace->elementToNodeMapping()->getElement(elementNo).nodeGlobalNo;
    std::vector<int> nodeGlobalNo;

    std::array<dof_no_t,FunctionSpaceType::nDofsPerElement()> dofsOfElement = functionSpace->getElementDofNosLocal(elementNo);
    for (typename std::array<dof_no_t,FunctionSpaceType::nDofsPerElement()>::const_iterator iter = dofsOfElement.begin(); iter != dofsOfElement.end(); iter++)
    {
      dof_no_t dofNo = *iter;
      if (dofNo % FunctionSpaceType::nDofsPerNode() == 0)
      {
        node_no_t nodeNo = dofNo / FunctionSpaceType::nDofsPerNode();
        nodeGlobalNo.push_back(nodeNo);
      }
    }

    // swap elements because VTK numbering for VTK_HEXAHEDRON (type 12) and VTK_QUAD (type 9) is different
    std::vector<int> permutatedNodeNos(nodeGlobalNo.begin(), nodeGlobalNo.end());

    if (permutatedNodeNos.size() >= 4)
    {
      /*int temp = permutatedNodeNos[2];
      permutatedNodeNos[2] = permutatedNodeNos[3];
      permutatedNodeNos[3] = temp;*/
      std::swap(permutatedNodeNos[2], permutatedNodeNos[3]);
    }

    if (permutatedNodeNos.size() >= 8)
    {
      /*int temp = permutatedNodeNos[6];
      permutatedNodeNos[6] = permutatedNodeNos[7];
      permutatedNodeNos[7] = temp;*/
      std::swap(permutatedNodeNos[6], permutatedNodeNos[7]);
    }

    connectivityValues.insert(connectivityValues.end(), permutatedNodeNos.begin(), permutatedNodeNos.end());
    /*for (int node : permutatedNodeNos)
    {
      connectivityValues.push_back(node);
    }*/
  }
  LOG(DEBUG) << ", now connectivityValues: " << connectivityValues;
}

}  // namespace OutputWriter
