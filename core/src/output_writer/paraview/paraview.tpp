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
                                          std::ofstream &file, bool binaryOutput, bool fixedFormat, bool onlyParallelDatasetElement)
{
 
  // here we have the type of the mesh with meshName (which is typedef to FunctionSpace)
  //typedef typename FieldVariableType::FunctionSpace FunctionSpace;

  // if only the "parallel dataset element" stub which is needed in the master files, should be written
  if (onlyParallelDatasetElement)
  {
    file << std::string(3, '\t') << "<PDataArray "
        << "Name=\"" << fieldVariable.name() << "\" "
        << "type=\"Float32\" "
        << "NumberOfComponents=\"" << fieldVariable.nComponents() << "\" ";

    if (binaryOutput)
    {
      file << "format=\"binary\" />" << std::endl;
    }
    else
    {
      file << "format=\"ascii\" />" << std::endl;
    }
  }
  else
  {
    // write normal data element
    file << std::string(4, '\t') << "<DataArray "
        << "Name=\"" << fieldVariable.name() << "\" "
        << "type=\"Float32\" "
        << "NumberOfComponents=\"" << fieldVariable.nComponents() << "\" ";

    const int nComponents = FieldVariableType::nComponents();
    std::string stringData;

    std::vector<double> values;
    std::array<std::vector<double>, nComponents> componentValues;

    // initialize the dofNosLocalNaturalOrdering vector of the meshPartition to be able to get the values in the natural ordering
    fieldVariable.functionSpace()->meshPartition()->initializeDofNosLocalNaturalOrdering(fieldVariable.functionSpace());

    // ensure that ghost values are in place
    fieldVariable.startGhostManipulation();

    // get all local values including ghosts for the components
    for (int componentNo = 0; componentNo < nComponents; componentNo++)
    {
      std::vector<double> retrievedLocalValues;
      fieldVariable.getValues(componentNo, fieldVariable.functionSpace()->meshPartition()->dofNosLocalNaturalOrdering(),
                              retrievedLocalValues);

      const int nDofsPerNode = FieldVariableType::FunctionSpace::nDofsPerNode();
      const node_no_t nNodesLocal = fieldVariable.functionSpace()->meshPartition()->nNodesLocalWithGhosts();

      // for Hermite only extract the non-derivative values
      componentValues[componentNo].resize(nNodesLocal);

      int index = 0;
      for (int i = 0; i < nNodesLocal; i++)
      {
        componentValues[componentNo][i] = retrievedLocalValues[index];
        index += nDofsPerNode;
      }
    }
    values.reserve(componentValues[0].size()*nComponents);

    // copy values in consecutive order (x y z x y z) to output
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
      file << "format=\"binary\" >" << std::endl;
    }
    else
    {
      stringData = Paraview::convertToAscii(values, fixedFormat);
      file << "format=\"ascii\" >" << std::endl;
    }
    file << std::string(5, '\t') << stringData << std::endl
      << std::string(4, '\t') << "</DataArray>" << std::endl;
  }
}
  
template<typename FieldVariableType>
void Paraview::writeParaviewPartitionFieldVariable(FieldVariableType &geometryField,
                                                   std::ofstream &file, bool binaryOutput, bool fixedFormat, bool onlyParallelDatasetElement)
{
  // if only the "parallel dataset element" stub which is needed in the master files, should be written
  if (onlyParallelDatasetElement)
  {
    file << std::string(3, '\t') << "<PDataArray "
        << "Name=\"partitioning\" "
        << "type=\"Int32\" "
        << "NumberOfComponents=\"1\" ";

    if (binaryOutput)
    {
      file << "format=\"binary\" />" << std::endl;
    }
    else
    {
      file << "format=\"ascii\" />" << std::endl;
    }
  }
  else
  {
    // write normal data element
    file << std::string(4, '\t') << "<DataArray "
        << "Name=\"partitioning\" "
        << "type=\"Int32\" "
        << "NumberOfComponents=\"1\" ";

    std::string stringData;

    // get own rank no
    int ownRankNoCommWorld;
    MPIUtility::handleReturnValue(MPI_Comm_rank(MPI_COMM_WORLD, &ownRankNoCommWorld));

    const node_no_t nNodesLocal = geometryField.functionSpace()->meshPartition()->nNodesLocalWithGhosts();

    std::vector<int> values(nNodesLocal, ownRankNoCommWorld);

    if (binaryOutput)
    {
      stringData = Paraview::encodeBase64(values);
      file << "format=\"binary\" >" << std::endl;
    }
    else
    {
      stringData = Paraview::convertToAscii(values, fixedFormat);
      file << "format=\"ascii\" >" << std::endl;
    }
    file << std::string(5, '\t') << stringData << std::endl
      << std::string(4, '\t') << "</DataArray>" << std::endl;
  }
}

};
