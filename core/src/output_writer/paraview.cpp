#include "output_writer/paraview.h"

#include <fstream>
#include <iostream>
#include <thread>
#include <chrono>

#include "easylogging++.h"
#include "base64.h"

#include <control/python_utility.h>
#include <control/petsc_utility.h>
#include <mesh/regular_fixed.h>
#include <mesh/rectilinear_fixed.h>
#include <mesh/nonrectilinear_fixed.h>
#include <mesh/deformable.h>
#include <mesh/mesh.h>

namespace OutputWriter
{

Paraview::Paraview(PyObject *settings) : Generic(settings)
{
}

void Paraview::writeSolution(Data::Data& data, int timeStepNo, double currentTime)
{
  if (!data.mesh())
  {
    LOG(FATAL) << "mesh is not set!";
  }
  const int dimension = data.mesh()->dimension();
  
  // solution and rhs vectors in mesh shape
  switch(dimension)
  {
  case 1:
    writeSolutionDim<1>(data, timeStepNo, currentTime);
    break;
  case 2:
    writeSolutionDim<2>(data, timeStepNo, currentTime);
    break;
  case 3:
    writeSolutionDim<3>(data, timeStepNo, currentTime);
    break;
  };
}

template <int dimension>
void Paraview::writeSolutionDim(Data::Data &data, int timeStepNo, double currentTime)
{
  LOG(TRACE) << "writeMesh<"<<dimension<<">()";
  
  if (std::dynamic_pointer_cast<Mesh::RegularFixed<dimension>>(data.mesh()) != NULL)
  {
    writeRectilinearGrid<Mesh::RegularFixed<dimension>>(data);
  }
  else if (std::dynamic_pointer_cast<Mesh::RectilinearFixed<dimension>>(data.mesh()) != NULL)
  {
    writeRectilinearGrid<Mesh::RectilinearFixed<dimension>>(data);
  }
  else if (std::dynamic_pointer_cast<Mesh::NonrectilinearFixed<dimension>>(data.mesh()) != NULL)
  {
    LOG(ERROR) << "not implemented";
  }
  else if (std::dynamic_pointer_cast<Mesh::Deformable<dimension>>(data.mesh()) != NULL)
  {
    LOG(ERROR) << "not implemented";
  }
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

template <class Mesh>
void Paraview::writeRectilinearGrid(Data::Data& data)
{
  // determine file name
  std::stringstream s;
  s<<filename_<<".vtr";
  std::string filename = s.str();

  // open file
  std::ofstream file(filename.c_str(), std::ios::out);
  
  if (!file.is_open())
  {
    // try to create directories
    if (filename.rfind("/") != std::string::npos)
    {
      // extract directory from filename
      std::string path = filename.substr(0, filename.rfind("/"));
      
      // create directory and wait until system has created it
      int ret = system((std::string("mkdir -p ")+path).c_str());
      if (ret != 0)
        LOG(WARNING) << "Creation of directory \""<<path<<"\" failed.";
      std::this_thread::sleep_for(std::chrono::milliseconds(500));
      
      file.clear();
      file.open(filename.c_str(), std::ios::out);
    }
  }
  
  if (!file.is_open())
  {
    LOG(WARNING) << "Could not open file \""<<filename<<"\" for writing!";
  }
  
  LOG(DEBUG) << "Write RectilinearGrid, file \""<<filename<<"\".";
  
  // determine values
  std::shared_ptr<Mesh> mesh = std::static_pointer_cast<Mesh>(data.mesh());
  
  // extent
  std::vector<int> extent = {0,0,0};   // number of nodes in x, y and z direction
  for (int i=0; i<mesh->dimension(); i++)
    extent[i] = mesh->nElements(i) + 1 - 1;   // nDOFS = nElements+1, value-1 because indices start with 0
  
  // coordinates of grid 
  std::vector<double> coordinates[3];
  int dimensionNo = 0;
  for (; dimensionNo<mesh->dimension(); dimensionNo++)
  {
    double meshWidth = mesh->meshWidth(dimensionNo);
    double nElements = mesh->nElements(dimensionNo);
    
    LOG(DEBUG) << "dimension "<<dimensionNo<<", meshWidth: "<<meshWidth<<", nElements: "<<nElements;
    
    coordinates[dimensionNo].resize(nElements+1);
    
    for(int nodeNo = 0; nodeNo < nElements+1; nodeNo++)
    {
      double coordinate = nodeNo * meshWidth;
      LOG(DEBUG) << "coordinate: "<<coordinate<<", nodeNo="<<nodeNo;
      coordinates[dimensionNo][nodeNo] = coordinate;
    }
  }
  
  // set other coordinates to 0
  for(; dimensionNo<3; dimensionNo++)
  {
    coordinates[dimensionNo].resize(1);
    coordinates[dimensionNo][0] = 0.0;
  }
  
  // name of value field
  std::string scalarsName = "Solution";
  
  bool binaryOutput = PythonUtility::getOptionBool(specificSettings_, "binaryOutput", true);
  bool fixedFormat = PythonUtility::getOptionBool(specificSettings_, "fixedFormat", true);
  
  // write file
  file << "<?xml version=\"1.0\"?>" << std::endl
    << "<VTKFile type=\"RectilinearGrid\" version=\"1.0\" byte_order=\"LittleEndian\">" << std::endl    // intel cpus are LittleEndian
    << std::string(1, '\t') << "<RectilinearGrid " 
      << "WholeExtent=\"" << "0 " << extent[0] << " 0 "<< extent[1] << " 0 " << extent[2] << "\"> " << std::endl     // dataset element
    << std::string(2, '\t') << "<Piece " 
      << "Extent=\"0 " << extent[0] << " 0 "<< extent[1] << " 0 " << extent[2] << "\"> " << std::endl
    << std::string(3, '\t') << "<PointData Scalars=\"" << scalarsName << "\">" << std::endl;
  if (binaryOutput)
  {
    file << std::string(4, '\t') << "<DataArray " 
        << "Name=\"" << scalarsName << "\" " 
        << "type=\"Float32\" " 
        << "NumberOfComponents=\"1\" "
        << "format=\"binary\" >" << std::endl
      << std::string(5, '\t') << encodeBase64(data.solution()) << std::endl
      << std::string(4, '\t') << "</DataArray>" << std::endl;
  }
  else
  {
    file << std::string(4, '\t') << "<DataArray " 
        << "Name=\"" << scalarsName << "\" " 
        << "type=\"Float32\" " 
        << "NumberOfComponents=\"1\" "
        << "format=\"ascii\" >" << std::endl
      << std::string(5, '\t') << convertToAscii(data.solution(), fixedFormat) << std::endl
      << std::string(4, '\t') << "</DataArray>" << std::endl;
  }
  
  file << std::string(3, '\t') << "</PointData>" << std::endl
    << std::string(3, '\t') << "<CellData>" << std::endl
    << std::string(3, '\t') << "</CellData>" << std::endl
    << std::string(3, '\t') << "<Coordinates>" << std::endl;
    
  if (binaryOutput)
  {
    file << std::string(4, '\t') << "<DataArray " 
        << "type=\"Float32\" " 
        << "NumberOfComponents=\"1\" "
        << "format=\"binary\" >" << std::endl
      << std::string(5, '\t') << encodeBase64(coordinates[0]) << std::endl
      << std::string(4, '\t') << "</DataArray>" << std::endl
      << std::string(4, '\t') << "<DataArray " 
        << "type=\"Float32\" " 
        << "NumberOfComponents=\"1\" "
        << "format=\"binary\" >" << std::endl
      << std::string(5, '\t') << encodeBase64(coordinates[1]) << std::endl
      << std::string(4, '\t') << "</DataArray>" << std::endl
      << std::string(4, '\t') << "<DataArray " 
        << "type=\"Float32\" " 
        << "NumberOfComponents=\"1\" "
        << "format=\"binary\" >" << std::endl
      << std::string(5, '\t') << encodeBase64(coordinates[2]) << std::endl
      << std::string(4, '\t') << "</DataArray>" << std::endl;
  }
  else
  {
    file << std::string(4, '\t') << "<DataArray " 
        << "type=\"Float32\" " 
        << "NumberOfComponents=\"1\" "
        << "format=\"ascii\" >" << std::endl
      << std::string(5, '\t') << convertToAscii(coordinates[0], fixedFormat) << std::endl
      << std::string(4, '\t') << "</DataArray>" << std::endl
      << std::string(4, '\t') << "<DataArray " 
        << "type=\"Float64\" " 
        << "NumberOfComponents=\"1\" "
        << "format=\"ascii\" >" << std::endl
      << std::string(5, '\t') << convertToAscii(coordinates[1], fixedFormat) << std::endl
      << std::string(4, '\t') << "</DataArray>" << std::endl
      << std::string(4, '\t') << "<DataArray " 
        << "type=\"Float64\" " 
        << "NumberOfComponents=\"1\" "
        << "format=\"ascii\" >" << std::endl
      << std::string(5, '\t') << convertToAscii(coordinates[2], fixedFormat) << std::endl
      << std::string(4, '\t') << "</DataArray>" << std::endl;
  }
  file << std::string(3, '\t') << "</Coordinates>" << std::endl
    << std::string(2, '\t') << "</Piece>" << std::endl
    << std::string(1, '\t') << "</RectilinearGrid>" << std::endl
    << "</VTKFile>"<<std::endl;
}


void Paraview::writeVTKMasterFile()
{
  /*
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
    
  os<<"</VTKFile>"<<endl;*/
}

void Paraview::writeVTKSlaveFile()
{

}


};