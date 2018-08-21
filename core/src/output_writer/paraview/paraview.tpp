#include "output_writer/paraview/paraview.h"

#include <fstream>
#include <iostream>
#include <thread>
#include <chrono>

#include "easylogging++.h"
#include "base64.h"

#include "output_writer/paraview/loop_output.h"
namespace OutputWriter
{

template<typename DataType>
void Paraview::write(DataType& data, int timeStepNo, double currentTime)
{
  // check if output should be written in this timestep and prepare filename
  if (!Generic::prepareWrite(data, timeStepNo, currentTime))
  {
    return;
  }

  // collect all available meshes
  std::set<std::string> meshNames;
  LoopOverTuple::loopCollectMeshNames<typename DataType::OutputFieldVariables>(data.getOutputFieldVariables(), meshNames);
  
  // loop over meshes and create a paraview file for each
  for (std::string meshName : meshNames)
  {
    // setup name of file
    std::stringstream filenameStart;
    if (meshNames.size() == 1)
      filenameStart << this->filename_;
    else
      filenameStart << this->filename_ << "_" << meshName;
   
    // loop over all field variables and output those that are associated with the mesh given by meshName
    ParaviewLoopOverTuple::loopOutput(data.getOutputFieldVariables(), data.getOutputFieldVariables(), meshName, filenameStart.str(), specificSettings_);
  }
}

template<typename FieldVariableType>
void Paraview::writeParaviewFieldVariable(FieldVariableType &fieldVariable, 
                                          std::ofstream &file, bool binaryOutput, bool fixedFormat)
{
 
  // here we have the type of the mesh with meshName (which is typedef to BasisOnMesh)
  typedef typename FieldVariableType::BasisOnMesh BasisOnMesh;

  file << std::string(4, '\t') << "<DataArray "
      << "Name=\"" << fieldVariable.name() << "\" "
      << "type=\"Float32\" "
      << "NumberOfComponents=\"" << fieldVariable.nComponents() << "\" ";
     
  const int nComponents = FieldVariableType::nComponents();
  std::string stringData;
  
  // for Hermite fields retrieve all data and interleave the components afterwards
  if (std::is_same<typename BasisOnMesh::BasisFunction, BasisFunction::Hermite>::value)
  {
    std::vector<double> values;
    std::array<std::vector<double>, nComponents> componentValues;
    
    for (int componentNo = 0; componentNo < nComponents; componentNo++)
    {
      fieldVariable.getValuesWithoutGhosts(componentNo, componentValues[componentNo], true);
    }
    values.reserve(componentValues[0].size()*nComponents);
    
    for (int i = 0; i < componentValues[0].size(); i++)
    {
      for (int componentNo = 0; componentNo < nComponents; componentNo++)
      {
        values.push_back(componentValues[componentNo][i]);
      }
    }
    
    if (binaryOutput)
    {
      stringData = Paraview::encodeBase64(values);
    }
    else 
    {
      stringData = Paraview::convertToAscii(values, fixedFormat);
    }
  }
  else    // not Hermite
  {
    // directly use the Petsc Vec
    if (binaryOutput)
    {
      stringData = Paraview::encodeBase64(fieldVariable.valuesLocal());  // wrong anyway
    }
    else 
    {
      stringData = Paraview::convertToAscii(fieldVariable.valuesLocal(), fixedFormat);
    }
  }
  
  if (binaryOutput)
  {
    file << "format=\"binary\" >" << std::endl;
  }
  else
  {
    file << "format=\"ascii\" >" << std::endl;
  }
  file << std::string(5, '\t') << stringData << std::endl
       << std::string(4, '\t') << "</DataArray>" << std::endl;
}
  
};