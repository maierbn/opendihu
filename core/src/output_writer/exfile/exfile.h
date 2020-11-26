#pragma once

#include <Python.h>  // has to be the first included header
#include <iostream>
#include <vector>

#include "control/types.h"
#include "output_writer/generic.h"

namespace OutputWriter
{

class Exfile : public Generic
{
public:

  //! constructor
  using Generic::Generic;

  //! write out solution to given filename, if timeStepNo is not -1, this value will be part of the filename
  template<typename DataType>
  void write(DataType &data, int timeStepNo = -1, double currentTime = -1, int callCountIncrement = 1);

private:

  //! output a cmgui "visualization.com" file that references the exfiles and can be used by Cmgui
  void outputComFile();

  struct FilenameWithElementAndNodeCount
  {
    std::string filename;  //< the filename 
    element_no_t nElements;  //< the number of elements 
    node_no_t nNodes;   //< the number of nodes
    std::string meshName;  //< the mesh name
    int dimensionality;   //< the dimensionality of the mesh, i.e. 1, 2 or 3
  };
  
  std::map<double,std::vector<FilenameWithElementAndNodeCount>> filenamesWithElementAndNodeCount_;   //< for a given simulation time the filenames without suffix of all previously output exelem files
};

} // namespace

#include "output_writer/exfile/exfile.tpp"
