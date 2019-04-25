#include "output_writer/paraview/poly_data_properties_for_mesh.h"

#include <Python.h>  // has to be the first included header
#include <iostream>
#include <vector>

#include "utility/vector_operators.h"
#include "control/types.h"

namespace OutputWriter
{

// output operator
std::ostream &operator<<(std::ostream &stream, PolyDataPropertiesForMesh rhs)
{
  stream << "(D: " << rhs.dimensionality
    << ", nPointsLocal: " << rhs.nPointsLocal << ", nCellsLocal: " << rhs.nCellsLocal
    << ", nPointsGlobal: " << rhs.nPointsGlobal << ", nCellsGlobal: " << rhs.nCellsGlobal
     << ", pointDataArrays:";
  for (auto pointDataArray : rhs.pointDataArrays)
  {
    stream << " [" << pointDataArray.first << "," << pointDataArray.second << "]";
  }
  stream << ")";
  return stream;
}

}  // namespace
