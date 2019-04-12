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
  // get connectivity information for the paraview vtu file
  // This method subdivides all mesh elements into cells with 2^D nodes (only for quadratic elements the cells are different from elements).
  // For every cell it adds the global node nos to the output vector connectivityValues
  int nNodesPerElement1D = FunctionSpace::FunctionSpaceBaseDim<1, BasisFunctionType>::nNodesPerElement();

  LOG(DEBUG) << "GetConnectivityValuesUnstructuredMesh, " << functionSpace->nElementsLocal() << " elements";

  // loop over elements of mesh
  for (element_no_t elementNo = 0; elementNo < functionSpace->nElementsLocal(); elementNo++)
  {
    // get global node nos of element from unstructured mesh data structure
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

    // extract the node nos for the cell
    if (D == 3)
    {
      // loop over cells of element (for quadratic 3D elements, there are 2x2x2=8 paraview cells per element)
      for (int k = 0; k < nNodesPerElement1D-1; k++)
      {
        for (int j = 0; j < nNodesPerElement1D-1; j++)
        {
          for (int i = 0; i < nNodesPerElement1D-1; i++)
          {
            std::vector<int> cellNodeGlobalNo(
            {
              nodeGlobalNo[k*nNodesPerElement1D*nNodesPerElement1D + j*nNodesPerElement1D + i],
              nodeGlobalNo[k*nNodesPerElement1D*nNodesPerElement1D + j*nNodesPerElement1D + (i+1)],
              nodeGlobalNo[k*nNodesPerElement1D*nNodesPerElement1D + (j+1)*nNodesPerElement1D + i],
              nodeGlobalNo[k*nNodesPerElement1D*nNodesPerElement1D + (j+1)*nNodesPerElement1D + (i+1)],
              nodeGlobalNo[(k+1)*nNodesPerElement1D*nNodesPerElement1D + j*nNodesPerElement1D + i],
              nodeGlobalNo[(k+1)*nNodesPerElement1D*nNodesPerElement1D + j*nNodesPerElement1D + (i+1)],
              nodeGlobalNo[(k+1)*nNodesPerElement1D*nNodesPerElement1D + (j+1)*nNodesPerElement1D + i],
              nodeGlobalNo[(k+1)*nNodesPerElement1D*nNodesPerElement1D + (j+1)*nNodesPerElement1D + (i+1)]
            });

            // swap elements because VTK numbering for VTK_HEXAHEDRON (type 12) is different
            std::vector<int> permutatedNodeNos(cellNodeGlobalNo.begin(), cellNodeGlobalNo.end());

            std::swap(permutatedNodeNos[2], permutatedNodeNos[3]);
            std::swap(permutatedNodeNos[6], permutatedNodeNos[7]);

            connectivityValues.insert(connectivityValues.end(), permutatedNodeNos.begin(), permutatedNodeNos.end());
          }
        }
      }
    }
    else if (D == 2)
    {
      // loop over cells of element (for quadratic 2D elements, there are 2x2=4 paraview cells per element)
      for (int j = 0; j < nNodesPerElement1D-1; j++)
      {
        for (int i = 0; i < nNodesPerElement1D-1; i++)
        {
          std::vector<int> cellNodeGlobalNo(
          {
            nodeGlobalNo[j*nNodesPerElement1D + i],
            nodeGlobalNo[j*nNodesPerElement1D + (i+1)],
            nodeGlobalNo[(j+1)*nNodesPerElement1D + i],
            nodeGlobalNo[(j+1)*nNodesPerElement1D + (i+1)]
          });

          // swap elements because VTK numbering for VTK_QUAD (type 9) is different
          std::vector<int> permutatedNodeNos(cellNodeGlobalNo.begin(), cellNodeGlobalNo.end());

          std::swap(permutatedNodeNos[2], permutatedNodeNos[3]);

          connectivityValues.insert(connectivityValues.end(), permutatedNodeNos.begin(), permutatedNodeNos.end());
        }
      }
    }
    else if (D == 1)
    {
      for (int i = 0; i < nNodesPerElement1D-1; i++)
      {
        std::vector<int> cellNodeGlobalNo(
        {
          nodeGlobalNo[i],
          nodeGlobalNo[i+1]
        });
        connectivityValues.insert(connectivityValues.end(), cellNodeGlobalNo.begin(), cellNodeGlobalNo.end());
      }
    }
  }
  LOG(DEBUG) << ", new connectivityValues: " << connectivityValues;
}

}  // namespace OutputWriter
