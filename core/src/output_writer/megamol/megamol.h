#pragma once

#include <Python.h>  // has to be the first included header
#include <iostream>
#include <vector>

#include "control/types.h"
#include "output_writer/generic.h"

namespace OutputWriter
{

struct BoundingBox
{
  Vec3 min, max;   // minimum,maximum value of bounding box
  bool initialized = false;   // if values are set
  BoundingBox();
};

class MegaMol : public Generic
{
public:

  //! constructor
  MegaMol(DihuContext context, PythonConfig specificSettings);

  //! write out solution to given filename, if timeStepNo is not -1, this value will be part of the filename
  template<typename DataType>
  void write(DataType &data, int timeStepNo = -1, double currentTime = -1);

private:

#if defined(HAVE_MEGAMOL) && defined(HAVE_ADIOS)
  void notifyMegaMol();
#endif

#ifdef HAVE_ADIOS
  std::shared_ptr<adios2::Engine> writer_;    //< adios writer
  std::shared_ptr<adios2::Variable<double>> boxVariable_ = nullptr; ///< variable that contains bounding box
  std::shared_ptr<adios2::Variable<int>> pCountVariable_ = nullptr; ///< variable that has number of particles/nodes
  std::shared_ptr<adios2::Variable<double>> globalRadiusVariable_ = nullptr; ///< variable that defines the radius of a sphere in megamol

  std::vector<double> boundingBoxValues_;   ///< bounding box
  int globalNumberOfNodes_;   ///< the "particle count" value for MegaMol
  std::string filenameSuffix_;   ///< the last part of the filename, either "_0" or "_1"
  std::string lastFilename_;   ///< the last filename that was written to
#endif

};

} // namespace

#include "output_writer/megamol/megamol.tpp"