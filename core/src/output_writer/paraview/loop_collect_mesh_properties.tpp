#include "output_writer/paraview/loop_collect_field_variables_names.h"

#include "output_writer/paraview/paraview_writer.h"
#include "output_writer/paraview/poly_data_properties_for_mesh.h"
#include "output_writer/paraview/get_connectivity_values_unstructured_mesh.h"

#include "function_space/00_function_space_base_dim.h"

#include <cstdlib>

namespace OutputWriter
{

namespace ParaviewLoopOverTuple
{
 
 /** Static recursive loop from 0 to number of entries in the tuple
 * Loop body
 */
template<typename FieldVariablesForOutputWriterType, int i>
inline typename std::enable_if<i < std::tuple_size<FieldVariablesForOutputWriterType>::value, void>::type
loopCollectMeshProperties(const FieldVariablesForOutputWriterType &fieldVariables, std::map<std::string,PolyDataPropertiesForMesh> &meshProperties
)
{
  LOG(DEBUG) << "loopCollectMeshProperties i=" << i << " type " << StringUtility::demangle(typeid(typename std::tuple_element<i,FieldVariablesForOutputWriterType>::type).name());

  // call what to do in the loop body
  if (collectMeshProperties<typename std::tuple_element<i,FieldVariablesForOutputWriterType>::type, FieldVariablesForOutputWriterType>(
        std::get<i>(fieldVariables), fieldVariables, meshProperties))
    return;
  
  // advance iteration to next tuple element
  loopCollectMeshProperties<FieldVariablesForOutputWriterType, i+1>(fieldVariables, meshProperties);
}
 
// current element is of pointer type (not vector)
template<typename CurrentFieldVariableType, typename FieldVariablesForOutputWriterType>
typename std::enable_if<!TypeUtility::isTuple<CurrentFieldVariableType>::value && !TypeUtility::isVector<CurrentFieldVariableType>::value, bool>::type
collectMeshProperties(CurrentFieldVariableType currentFieldVariable, const FieldVariablesForOutputWriterType &fieldVariables,
                           std::map<std::string,PolyDataPropertiesForMesh> &meshProperties)
{
  assert(currentFieldVariable != nullptr);
  assert(currentFieldVariable->functionSpace());
  std::string meshName = currentFieldVariable->functionSpace()->meshName();
  LOG(DEBUG) << "field variable \"" << currentFieldVariable->name() << "\", mesh \"" << meshName << "\".";

  /*
  int dimensionality;    ///< D=1: object is a VTK "Line", D=2, D=3: object should be represented by an unstructured grid
  global_no_t nPoints;   ///< the number of points needed for representing the mesh
  global_no_t nCells;    ///< the number of VTK "cells", i.e. "Lines" or "Polys", which is the opendihu number of "elements"
  std::vector<node_no_t> nNodesLocalWithGhosts;   ///< local number of nodes including ghosts, for all dimensions

  std::vector<std::pair<std::string,int>> pointDataArrays;   ///< <name,nComponents> of PointData DataArray elements
  */
  int dimensionality = currentFieldVariable->functionSpace()->dim();
  meshProperties[meshName].dimensionality = dimensionality;
  meshProperties[meshName].nPointsLocal = currentFieldVariable->functionSpace()->nNodesLocalWithGhosts();
  meshProperties[meshName].nPointsGlobal = currentFieldVariable->functionSpace()->nNodesGlobal();

  typedef typename CurrentFieldVariableType::element_type::FunctionSpace FunctionSpaceType;
  typedef typename FunctionSpaceType::BasisFunction BasisFunction;

  int nNodesPerElement1D = FunctionSpace::FunctionSpaceBaseDim<1, BasisFunction>::nNodesPerElement();
  int nCellsPerElement = pow(nNodesPerElement1D-1,dimensionality);

  meshProperties[meshName].nCellsGlobal = currentFieldVariable->functionSpace()->nElementsGlobal() * nCellsPerElement;
  meshProperties[meshName].nCellsLocal = currentFieldVariable->functionSpace()->nElementsLocal() * nCellsPerElement;

  meshProperties[meshName].nNodesLocalWithGhosts.resize(dimensionality);

  for (int dimensionIndex = 0; dimensionIndex < dimensionality; dimensionIndex++)
  {
    meshProperties[meshName].nNodesLocalWithGhosts[dimensionIndex] = currentFieldVariable->functionSpace()->meshPartition()->nNodesLocalWithGhosts(dimensionIndex);
  }

  // for unstructured grids add node nos per element for connectivity array
  if (meshProperties[meshName].unstructuredMeshConnectivityValues.empty())
  {
    GetConnectivityValuesUnstructuredMesh<typename CurrentFieldVariableType::element_type::FunctionSpace>::get(currentFieldVariable->functionSpace(), meshProperties[meshName].unstructuredMeshConnectivityValues);
  }


  if (!currentFieldVariable->isGeometryField())
  {
    meshProperties[meshName].pointDataArrays.push_back(std::pair<std::string,int>(currentFieldVariable->name(),currentFieldVariable->nComponents()));
  }

  return false;  // do not break iteration over field variables
}

// element i is of vector type
template<typename VectorType, typename FieldVariablesForOutputWriterType>
typename std::enable_if<TypeUtility::isVector<VectorType>::value, bool>::type
collectMeshProperties(VectorType currentFieldVariableVector, const FieldVariablesForOutputWriterType &fieldVariables,
                           std::map<std::string,PolyDataPropertiesForMesh> &meshProperties)
{
  for (auto& currentFieldVariable : currentFieldVariableVector)
  {
    // call function on all vector entries
    if (collectMeshProperties<typename VectorType::value_type,FieldVariablesForOutputWriterType>(currentFieldVariable, fieldVariables, meshProperties))
      return true; // break iteration
  }
  return false;  // do not break iteration 
}

// element i is of tuple type
template<typename TupleType, typename FieldVariablesForOutputWriterType>
typename std::enable_if<TypeUtility::isTuple<TupleType>::value, bool>::type
collectMeshProperties(TupleType currentFieldVariableTuple, const FieldVariablesForOutputWriterType &fieldVariables,
                           std::map<std::string,PolyDataPropertiesForMesh> &meshProperties)
{
  // call for tuple element
  loopCollectMeshProperties<TupleType>(currentFieldVariableTuple, meshProperties);
 
  return false;  // do not break iteration 
}

}  // namespace ParaviewLoopOverTuple
}  // namespace OutputWriter
