#include "output_writer/paraview/loop_collect_field_variables_names.h"

#include "output_writer/paraview/paraview_writer.h"
#include "output_writer/paraview/poly_data_properties_for_mesh.h"
#include "output_writer/paraview/get_connectivity_values_unstructured_mesh.h"

#include "function_space/00_function_space_base_dim.h"

#include <cstdlib>
#include "field_variable/field_variable.h"

namespace OutputWriter
{

namespace ParaviewLoopOverTuple
{
 
 /** Static recursive loop from 0 to number of entries in the tuple
 * Loop body
 */
template<typename FieldVariablesForOutputWriterType, int i>
inline typename std::enable_if<i < std::tuple_size<FieldVariablesForOutputWriterType>::value, void>::type
loopCollectMeshProperties(const FieldVariablesForOutputWriterType &fieldVariables, std::map<std::string,PolyDataPropertiesForMesh> &meshProperties,
                          std::vector<std::string> &meshNamesVector
)
{
  //LOG(DEBUG) << "loopCollectMeshProperties i=" << i << " type " << StringUtility::demangle(typeid(typename std::tuple_element<i,FieldVariablesForOutputWriterType>::type).name());

  // call what to do in the loop body
  if (collectMeshProperties<typename std::tuple_element<i,FieldVariablesForOutputWriterType>::type, FieldVariablesForOutputWriterType>(
        std::get<i>(fieldVariables), fieldVariables, meshProperties, meshNamesVector, i))
    return;
  
  // advance iteration to next tuple element
  loopCollectMeshProperties<FieldVariablesForOutputWriterType, i+1>(fieldVariables, meshProperties, meshNamesVector);
}
 
// current element is of pointer type (not vector)
template<typename CurrentFieldVariableType, typename FieldVariablesForOutputWriterType>
typename std::enable_if<!TypeUtility::isTuple<CurrentFieldVariableType>::value && !TypeUtility::isVector<CurrentFieldVariableType>::value && !Mesh::isComposite<CurrentFieldVariableType>::value, bool>::type
collectMeshProperties(CurrentFieldVariableType currentFieldVariable, const FieldVariablesForOutputWriterType &fieldVariables,
                           std::map<std::string,PolyDataPropertiesForMesh> &meshProperties, std::vector<std::string> &meshNamesVector, int i)
{
  if (currentFieldVariable == nullptr)
  {
    LOG(FATAL) << "In collectMeshProperties, currentFieldVariable is nullptr.\n"
      << " meshProperties: " << meshProperties
      << ", fielVariableType: " << StringUtility::demangle(typeid(CurrentFieldVariableType).name())
      << " fieldVariables: " << fieldVariables << ", i: " << i;
  }
  assert(currentFieldVariable != nullptr);
  if (!currentFieldVariable->functionSpace())
  {
    LOG(DEBUG) << "In collectMeshProperties, currentFieldVariable->functionSpace() is nullptr.\n"
      << " meshProperties: " << meshProperties
      << ", fielVariableType: " << StringUtility::demangle(typeid(CurrentFieldVariableType).name())
      << " fieldVariables: " << fieldVariables << ", i: " << i;
    return false;
  }
  assert(currentFieldVariable->functionSpace());
  std::string meshName = currentFieldVariable->functionSpace()->meshName();

  // if meshName is not already in meshNamesVector, add it there
  std::set<std::string> meshNamesSet(meshNamesVector.begin(),meshNamesVector.end());
  if (meshNamesSet.find(meshName) == meshNamesSet.end())
    meshNamesVector.push_back(meshName);
  //LOG(DEBUG) << "field variable \"" << currentFieldVariable->name() << "\", mesh \"" << meshName << "\".";

  /*
  int dimensionality;    ///< D=1: object is a VTK "Line", D=2, D=3: object should be represented by an unstructured grid
  global_no_t nPoints;   ///< the number of points needed for representing the mesh
  global_no_t nCells;    ///< the number of VTK "cells", i.e. "Lines" or "Polys", which is the opendihu number of "elements"
  std::vector<node_no_t> nNodesLocalWithGhosts;   ///< local number of nodes including ghosts, for all dimensions

  std::vector<PointData> pointDataArrays;   ///< <name,nComponents> of PointData DataArray elements
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
    PolyDataPropertiesForMesh::DataArrayName dataArrayName;
    dataArrayName.name = currentFieldVariable->name();
    dataArrayName.nComponents = currentFieldVariable->nComponents();
    dataArrayName.componentNames = std::vector<std::string>(currentFieldVariable->componentNames().begin(), currentFieldVariable->componentNames().end());

    meshProperties[meshName].pointDataArrays.push_back(dataArrayName);
  }

  return false;  // do not break iteration over field variables
}

// element i is of vector type
template<typename VectorType, typename FieldVariablesForOutputWriterType>
typename std::enable_if<TypeUtility::isVector<VectorType>::value, bool>::type
collectMeshProperties(VectorType currentFieldVariableGradient, const FieldVariablesForOutputWriterType &fieldVariables,
                      std::map<std::string,PolyDataPropertiesForMesh> &meshProperties,
                      std::vector<std::string> &meshNamesVector, int i)
{
  for (auto& currentFieldVariable : currentFieldVariableGradient)
  {
    // call function on all vector entries
    if (collectMeshProperties<typename VectorType::value_type,FieldVariablesForOutputWriterType>(
          currentFieldVariable, fieldVariables, meshProperties, meshNamesVector, i))
      return true; // break iteration
  }
  return false;  // do not break iteration 
}

// element i is of tuple type
template<typename TupleType, typename FieldVariablesForOutputWriterType>
typename std::enable_if<TypeUtility::isTuple<TupleType>::value, bool>::type
collectMeshProperties(TupleType currentFieldVariableTuple, const FieldVariablesForOutputWriterType &fieldVariables,
                      std::map<std::string,PolyDataPropertiesForMesh> &meshProperties,
                      std::vector<std::string> &meshNamesVector, int i)
{
  // call for tuple element
  loopCollectMeshProperties<TupleType>(currentFieldVariableTuple, meshProperties, meshNamesVector);
 
  return false;  // do not break iteration 
}

// element i is a field variables with Mesh::CompositeOfDimension<D>
template<typename CurrentFieldVariableType, typename FieldVariablesForOutputWriterType>
typename std::enable_if<Mesh::isComposite<CurrentFieldVariableType>::value, bool>::type
collectMeshProperties(CurrentFieldVariableType currentFieldVariable, const FieldVariablesForOutputWriterType &fieldVariables,
                      std::map<std::string,PolyDataPropertiesForMesh> &meshProperties,
                      std::vector<std::string> &meshNamesVector, int i)
{
  const int D = CurrentFieldVariableType::element_type::FunctionSpace::dim();
  typedef typename CurrentFieldVariableType::element_type::FunctionSpace::BasisFunction BasisFunctionType;
  typedef FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<D>, BasisFunctionType> SubFunctionSpaceType;
  const int nComponents = CurrentFieldVariableType::element_type::nComponents();

  typedef FieldVariable::FieldVariable<SubFunctionSpaceType, nComponents> SubFieldVariableType;

  std::vector<std::shared_ptr<SubFieldVariableType>> subFieldVariables;
  currentFieldVariable->getSubFieldVariables(subFieldVariables);

  for (auto& currentSubFieldVariable : subFieldVariables)
  {
    // call function on all vector entries
    if (collectMeshProperties<std::shared_ptr<SubFieldVariableType>,FieldVariablesForOutputWriterType>(
          currentSubFieldVariable, fieldVariables, meshProperties, meshNamesVector, i))
      return true;
  }

  return false;  // do not break iteration
}
}  // namespace ParaviewLoopOverTuple
}  // namespace OutputWriter
