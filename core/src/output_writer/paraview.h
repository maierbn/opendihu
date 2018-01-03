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
  Paraview(PyObject *specificSettings);
 
private:
 
  //! write out solution to given filename
  void writeSolution(Data::Data &data, int timeStepNo, double currentTime);
 
  //! write out solution templated by dimension 
  template <int dimension>
  void writeSolutionDim(Data::Data &data, int timeStepNo, double currentTime);
  
  //! write serial vtkRectilinearGrid file (structured, suffix *.vtr)
  template <class Mesh>
  void writeRectilinearGrid(Data::Data& data);
 
  //! write serial vtkStructuredGrid file (structured, suffix *.vts)
  template <class Mesh>
  void writeStructuredGrid(Data::Data& data);
  
  //! open file given by filename, create directory if necessary
  std::ofstream openFile(std::string filename);
 
  //! create *.pvt VTK master file that is a header file for multiple
  /** This writes the master file of the parallel output, should only be done
    * by one processor. */
  void writeVTKMasterFile();
 
  void writeVTKSlaveFile();
  
  //! encode a Petsc vector in Base64
  std::string encodeBase64(Vec &vector);
  
  //! encode a std::vector as base64
  std::string encodeBase64(std::vector<double> &vector);
  
  //! convert to a string with space separated values
  std::string convertToAscii(Vec &vector, bool humanReadable);
  std::string convertToAscii(std::vector<double> &vector, bool humanReadable);
};

};