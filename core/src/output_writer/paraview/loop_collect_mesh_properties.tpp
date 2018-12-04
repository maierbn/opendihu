#include "output_writer/paraview/loop_collect_field_variables_names.h"

#include "output_writer/paraview/paraview_writer.h"
#include "output_writer/paraview/poly_data_properties_for_mesh.h"

#include <cstdlib>

namespace OutputWriter
{

namespace ParaviewLoopOverTuple
{
 
 /** Static recursive loop from 0 to number of entries in the tuple
 * Loop body
 */
template<typename OutputFieldVariablesType, int i>
inline typename std::enable_if<i < std::tuple_size<OutputFieldVariablesType>::value, void>::type
loopCollectMeshProperties(const OutputFieldVariablesType &fieldVariables, std::map<std::string,PolyDataPropertiesForMesh> &meshProperties
)
{
  // call what to do in the loop body
  if (collectMeshProperties<typename std::tuple_element<i,OutputFieldVariablesType>::type, OutputFieldVariablesType>(
        std::get<i>(fieldVariables), fieldVariables, meshProperties))
    return;
  
  // advance iteration to next tuple element
  loopCollectMeshProperties<OutputFieldVariablesType, i+1>(fieldVariables, meshProperties);
}
 
// current element is of pointer type (not vector)
template<typename CurrentFieldVariableType, typename OutputFieldVariablesType>
typename std::enable_if<!TypeUtility::isTuple<CurrentFieldVariableType>::value && !TypeUtility::isVector<CurrentFieldVariableType>::value, bool>::type
collectMeshProperties(CurrentFieldVariableType currentFieldVariable, const OutputFieldVariablesType &fieldVariables,
                           std::map<std::string,PolyDataPropertiesForMesh> &meshProperties)
{
  std::string meshName = currentFieldVariable->functionSpace()->meshName();


  /*
  int dimensionality;    ///< D=1: object is a VTK "Line", D=2, D=3: object is a VTK "Poly"
  global_no_t nPoints;   ///< the number of points needed for representing the mesh
  global_no_t nCells;    ///< the number of VTK "cells", i.e. "Lines" or "Polys", which is the opendihu number of "elements"

  std::vector<std::pair<std::string,int>> pointDataArrays;   ///< <name,nComponents> of PointData DataArray elements
  */
  meshProperties[meshName].dimensionality = currentFieldVariable->functionSpace()->dim();
  meshProperties[meshName].nPointsLocal = currentFieldVariable->functionSpace()->nNodesLocalWithGhosts();
  meshProperties[meshName].nCellsLocal = currentFieldVariable->functionSpace()->nElementsLocal();
  meshProperties[meshName].nPointsGlobal = currentFieldVariable->functionSpace()->nNodesGlobal();
  meshProperties[meshName].nCellsGlobal = currentFieldVariable->functionSpace()->nElementsGlobal();

  if (!currentFieldVariable->isGeometryField())
  {
    meshProperties[meshName].pointDataArrays.push_back(std::pair<std::string,int>(currentFieldVariable->name(),currentFieldVariable->nComponents()));
  }

  return false;  // do not break iteration over field variables
}

// element i is of vector type
template<typename VectorType, typename OutputFieldVariablesType>
typename std::enable_if<TypeUtility::isVector<VectorType>::value, bool>::type
collectMeshProperties(VectorType currentFieldVariableVector, const OutputFieldVariablesType &fieldVariables,
                           std::map<std::string,PolyDataPropertiesForMesh> &meshProperties)
{
  for (auto& currentFieldVariable : currentFieldVariableVector)
  {
    // call function on all vector entries
    if (collectMeshProperties<typename VectorType::value_type,OutputFieldVariablesType>(currentFieldVariable, fieldVariables, meshProperties))
      return true; // break iteration
  }
  return false;  // do not break iteration 
}

// element i is of tuple type
template<typename TupleType, typename OutputFieldVariablesType>
typename std::enable_if<TypeUtility::isTuple<TupleType>::value, bool>::type
collectMeshProperties(TupleType currentFieldVariableTuple, const OutputFieldVariablesType &fieldVariables,
                           std::map<std::string,PolyDataPropertiesForMesh> &meshProperties)
{
  // call for tuple element
  loopCollectMeshProperties<TupleType>(currentFieldVariableTuple, meshProperties);
 
  return false;  // do not break iteration 
}

};  //namespace ParaviewLoopOverTuple
};  //namespace OutputWriter
