#pragma once

#include <iostream>
#include <vector>

#include "control/types.h"
#include "output_writer/generic.h"

namespace OutputWriter
{   

class Paraview : public Generic
{
public:
 
  //! constructor
  Paraview(PyObject *settings);
 
private:
 
  //! write serial vtkRectilinearGrid file (structured, suffix *.vts)
  template <class Mesh>
  void writeRectilinearGrid(Data::Data& data);
 
  //! create *.pvt VTK master file that is a header file for multiple
  /** This writes the master file of the parallel output, should only be done
    * by one processor. */
  void writeVTKMasterFile();
 
  void writeVTKSlaveFile();
  
  void writeSolution(Data::Data &data);
  
  template <int dimension>
  void writeSolutionDim(Data::Data &data);
  
  //! encode a Petsc vector in Base64
  std::string encodeBase64(Vec &vector);
  
  //! encode a std::vector as base64
  std::string encodeBase64(std::vector<double> &vector);
  
  //! convert to a string with space separated values
  std::string convertToAscii(Vec &vector, bool humanReadable);
  std::string convertToAscii(std::vector<double> &vector, bool humanReadable);
  
  PyObject *settings_;
};

};