#pragma once

#include <Python.h>  // has to be the first included header
#include <iostream>
#include <vector>

#include "control/types.h"
#include "output_writer/generic.h"

namespace OutputWriter
{

/**
 * The properties of a VTK PolyData object for one opendihu mesh.
 * This information on all ranks is sufficient to collectively compute file pointers of the output file for MPI file I/O.
 */
struct PolyDataPropertiesForMesh
{
  int dimensionality;    ///< D=1: object is a VTK "Line", D=2, D=3: object is a VTK "Poly"
  global_no_t nPoints;   ///< the number of points needed for representing the mesh
  global_no_t nCells;    ///< the number of VTK "cells", i.e. "Lines" or "Polys", which is the opendihu number of "elements"

  std::vector<std::pair<std::string,int>> pointDataArrays;   ///< <name,nComponents> of PointData DataArray elements
};

};
