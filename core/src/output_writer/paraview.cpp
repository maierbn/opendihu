#include "output_writer/paraview.h"

#include <fstream>
#include <iostream>
#include <thread>
#include <chrono>

#include "easylogging++.h"
#include "base64.h"

#include <utility/python_utility.h>
#include <utility/petsc_utility.h>
#include <mesh/regular_fixed.h>
#include <mesh/structured_deformable.h>
#include <mesh/unstructured_deformable.h>
#include <mesh/mesh.h>

namespace OutputWriter
{

Paraview::Paraview(PyObject *settings) : Generic(settings)
{
}

std::string Paraview::encodeBase64(Vec &vector)
{
  int vectorSize = 0;
  VecGetSize(vector, &vectorSize);
  
  std::vector<int> indices(vectorSize);
  std::iota (indices.begin(), indices.end(), 0);    // fill with increasing numbers: 0,1,2,...
  std::vector<double> values(vectorSize);
  VecGetValues(vector, vectorSize, indices.data(), values.data());
  
  return encodeBase64(values);
}

std::string Paraview::convertToAscii(Vec &vector, bool fixedFormat)
{
  std::vector<double> vectorValues;
  PetscUtility::getVectorEntries(vector, vectorValues);
  
  return convertToAscii(vectorValues, fixedFormat);
}

std::string Paraview::encodeBase64(std::vector<double> &vector)
{
  int rawLength = vector.size()*sizeof(float);
  int encodedLength = Base64::EncodedLength(rawLength);
  
  char raw[rawLength];
  for (unsigned int i=0; i<vector.size(); i++)
  {
    union {
      float d;
      char c[4];
    };
    d = vector[i];
    memcpy(raw+i*sizeof(float), c, 4);
  }
  
  char encoded[encodedLength+1];
  //Base64::Encode(reinterpret_cast<char *>(vector.data()), rawLength, encoded, encodedLength);
  bool success = Base64::Encode(raw, rawLength, encoded, encodedLength);
  if (!success)
    LOG(WARNING) << "encoding failed";
  
  encoded[encodedLength] = '\0';
  
  return std::string(encoded);
}

std::string Paraview::convertToAscii(std::vector<double> &vector, bool fixedFormat)
{
  std::stringstream result;
  for(auto value : vector)
  {
    if(fixedFormat)
    {
      result << std::setw(16) << std::scientific << value << " ";
    }
    else
    {
      result << value << " ";
    }
  }
  return result.str();
}
/*
void Paraview::writeVTKMasterFile()
{
  
  file<<"<?xml version=\"1.0\"?>"<<endl
    <<"<VTKFile type=\"PRectilinearGrid\">"<<endl
    <<"<PRectilinearGrid WholeExtent=\"0 "<<comm->get_global_dimension()[0]<<" 0 "
      <<comm->get_global_dimension()[1]<<" 0 0\" GhostLevel=\"0\">"
      <<endl
    <<"<PCoordinates>"<<endl
    <<"<PDataArray type=\"Float64\"/>"<<endl
    <<"<PDataArray type=\"Float64\"/>"<<endl
    <<"<PDataArray type=\"Float64\"/>"<<endl
    <<"</PCoordinates>"<<endl;

    MultiIndexType n_subareas = comm->get_n_subareas();

    int x_begin = 0;
    int y_begin = 0;
    
    MultiIndexType local_dimension;
    
    for(int y = 0; y < n_subareas[1]; y++)      // iterate over row (y-direction)
    {
        x_begin = 0;
        for(int x = 0; x < n_subareas[0]; x++)      // iterate over column (x-direction)
        {
            local_dimension = comm->get_local_dimension(x, y);
            
            os<<"<Piece Extent=\""
              <<x_begin<<" "<<x_begin+local_dimension[0]<<" "
              <<y_begin<<" "<<y_begin+local_dimension[1]<<" 0 0\" "
              <<"Source=\"field_"<<step<<"_processor_"<<y<<"_"<<x<<".vtr\"/>"<<endl;
              
            x_begin += local_dimension[0];
        }
        y_begin += local_dimension[1];
    }

  os<<"<PPointData Vectors=\"field\" Scalars=\"p, vorticity, stream\">"<<endl;
  os<<"<PDataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"field\"     format=\"ascii\"/>"<<endl;
  os<<"<PDataArray type=\"Float64\"                          Name=\"p\"         format=\"ascii\"/>"<<endl;
  os<<"<PDataArray type=\"Float64\"                          Name=\"vorticity\" format=\"ascii\"/>"<<endl;
  os<<"<PDataArray type=\"Float64\"                          Name=\"stream\"    format=\"ascii\"/>"<<endl;
    os<<"<PDataArray type=\"Float64\"                          Name=\"rhs\"       format=\"ascii\"/>"<<endl;
    os<<"<PDataArray type=\"Float64\"                          Name=\"partitioning\"    format=\"ascii\"/>"<<endl;
  os<<"</PPointData>"<<endl;

  os<<"</PRectilinearGrid>"<<endl;
    
    os<<"<PUnstructuredGrid GhostLevel=\"0\">"<<endl;
    os<<"<PPointData></PPointData>"<<endl;
    os<<"<PCellData>"<<endl;
    os<<"<DataArray type=\"Float32\" Name=\"particles\" format=\"ascii\">"<<endl;
    os<<"</DataArray>"<<endl;
    os<<"</PCellData>"<<endl;
    os<<"<PPoints>"<<endl;
    os<<"<PDataArray NumberOfComponents=\"3\">"<<endl;
    os<<"</PDataArray>"<<endl;
    os<<"</PPoints>"<<endl;
    os<<"<Piece Source=\"particles_"<<step<<".vtu\"/>"<<endl;

    os<<"</PUnstructuredGrid> "<<endl;
    
  os<<"</VTKFile>"<<endl;
}

void Paraview::writeVTKSlaveFile()
{

}
*/

};