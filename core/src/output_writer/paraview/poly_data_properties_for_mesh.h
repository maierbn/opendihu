#pragma once

#include <Python.h>  // has to be the first included header
#include <iostream>
#include <vector>

#include "control/types.h"

namespace OutputWriter
{

/**
 * The properties of a VTK PolyData object for one opendihu mesh.
 * This information on all ranks is sufficient to collectively compute file pointers of the output file for MPI file I/O.
 */
struct PolyDataPropertiesForMesh
{
  int dimensionality;           //< D=1: object is a VTK "Line", D=2, D=3: object should be represented by an unstructured grid
  global_no_t nPointsLocal;     //< the number of points needed for representing the mesh, local value of rank
  global_no_t nCellsLocal;      //< the number of VTK "cells", i.e. "Lines" or "Polys", which is the opendihu number of "elements", local value of rank
  global_no_t nPointsGlobal;    //< the number of points needed for representing the mesh, global value of all rank
  global_no_t nCellsGlobal;     //< the number of VTK "cells", i.e. "Lines" or "Polys", which is the opendihu number of "elements", global value of all ranks
  std::vector<node_no_t> nNodesLocalWithGhosts;         //< local number of nodes including ghosts, for all dimensions
  std::vector<int> unstructuredMeshConnectivityValues;  //< the values that will be written as "connectivity" array in the vtu file, when outputting UnstructuredDeformable meshes, these values are set directly

  struct DataArrayName
  {
    std::string name;                         //< name of the field variable
    int nComponents;                          //< number of components
    std::vector<std::string> componentNames;  //< names of the components
  };
  std::vector<PolyDataPropertiesForMesh::DataArrayName> pointDataArrays; //< names of PointData DataArray elements
};

// output operator
std::ostream &operator<<(std::ostream &stream, PolyDataPropertiesForMesh rhs);

}  // namespace
