//#include <vtkXMLUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkHexahedron.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <cmath>

#include "easylogging++.h"

//#include <experimental/filesystem>
//using namespace std::experimental::filesystem;

#include "boost/filesystem.hpp"
#include <iostream>
using namespace boost::filesystem;

// the data for one timestep
struct EmgData
{
  int width;
  int height;
  std::vector<double> phi;
};

// the data of one directory, i.e. for multiple timesteps
struct SimulationData
{
  std::vector<double> times;       //< time points of the parsed values
  std::vector<EmgData> values;   //< the values for a single time step, each, corresponding to times
};

void parseData(std::string filename, EmgData &data, double &t)
{
  // Create a reader
  vtkSmartPointer<vtkXMLUnstructuredGridReader> reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
 
  reader->SetFileName(filename.c_str());
  reader->Update();
  
  if (reader->GetNumberOfPointArrays() == 0)
  {
    LOG(ERROR) << "File \"" << filename << "\" contains no data.";
  }

  vtkSmartPointer<vtkDataSet> dataSet = reader->GetOutputAsDataSet(0);
  vtkSmartPointer<vtkFieldData> fieldData = dataSet->GetFieldData();
  
  /*
  std::cout << "file \"" << filename << "\" contains " << reader->GetNumberOfPointArrays() << " point arrays, "
    << reader->GetNumberOfCellArrays() << "," << reader->GetNumberOfColumnArrays() << std::endl;
  std::cout << " field has " << fieldData->GetNumberOfArrays() << " arrays.";
  */
  
  // get stored time
  if (fieldData->GetNumberOfArrays() > 0)
  {
    vtkSmartPointer<vtkDataArray> array = fieldData->GetArray(0);
    double *value = array->GetTuple(0);
    t = *value;
  }
  else 
  {
    LOG(WARNING) << "File \"" << filename << "\" contains no time.";
  }
  /*
  // get stored velocity values
  vtkSmartPointer<vtkPointData> pointData = dataSet->GetPointData();
  for (int arrayNo = 0; arrayNo < pointData->GetNumberOfArrays(); arrayNo++)
  {
    vtkSmartPointer<vtkDoubleArray> array = vtkDoubleArray::SafeDownCast(pointData->GetAbstractArray(arrayNo));
    int nEntries = array->GetNumberOfTuples();
    //int nComponents = array->GetNumberOfComponents();
    //std::cout << "  " << array->GetName() << ", nEntries: " << nEntries << ", nComponents: " << nComponents << std::endl;
    
    if (std::string(array->GetName()) != "velocity")
      continue;
    
    meshData.u.resize(nEntries);
    meshData.v.resize(nEntries);
    
    // loop over values
    for (int entryNo = 0; entryNo < nEntries; entryNo++)
    {
      std::array<double,3> values;
      array->GetTuple(entryNo, values.data());
      meshData.u[entryNo] = values[0];
      meshData.v[entryNo] = values[1];
    }
  }*/
}

void parseDirectory(std::string directory, SimulationData &simulationData)
{
  // get files in directory
  std::vector<std::string> filenames;
  
  directory_iterator end_iter; // default construction yields past-the-end
  for (directory_iterator iter(directory); iter != end_iter; ++iter)
  {
    if (is_directory(iter->status()))
      continue;
    
    if (is_regular_file(iter->status()))
    {
      std::string filename(iter->path().string());
      
      if (filename.find(".vtu") != std::string::npos)
      {
        filenames.push_back(filename);
      }
    }
  }
  
  // sort files by filename
  std::sort(filenames.begin(), filenames.end());
      
  //std::cout << "Parse data in \"" << directory << "\"." << std::endl;
  
  // loop over filenames and parse files
  for (std::vector<std::string>::iterator iter = filenames.begin(); iter != filenames.end(); iter++)
  {
    std::string filename = *iter;
    //std::cout << "  File " << filename << std::endl;
    double t;
    simulationData.values.emplace_back();
    parseData(filename, simulationData.values.back(), t);
    simulationData.times.push_back(t);
  }
}

int main(int argc, char *argv[])
{
  // parse command line arguments
  if (argc != 2)
  {
    std::cout << "usage: " << argv[0] << " <directory with data>" << std::endl;
    exit(-1);
  }
  
  std::string directory = argv[1];
  
  // parse data in the directory
  EmgData data;
  parseDirectory(directory, data);
  
  return 0;
}
       
