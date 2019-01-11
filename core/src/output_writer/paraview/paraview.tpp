#include "output_writer/paraview/paraview.h"

#include <fstream>
#include <iostream>
#include <thread>
#include <chrono>

#include "easylogging++.h"
#include "base64.h"

#include "output_writer/paraview/loop_output.h"
#include "output_writer/paraview/loop_collect_mesh_properties.h"
#include "output_writer/paraview/poly_data_properties_for_mesh.h"

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

  std::set<std::string> combined1DMeshes;
  if (combineFiles_)
  {
    // create a PolyData file that combines all 1D meshes into one file
    writePolyDataFile<typename DataType::OutputFieldVariables>(data.getOutputFieldVariables(), combined1DMeshes);
  }

  // output normal files, parallel or if combineFiles_, only the 2D and 3D meshes, combined

  // collect all available meshes
  std::set<std::string> meshNames;
  LoopOverTuple::loopCollectMeshNames<typename DataType::OutputFieldVariables>(data.getOutputFieldVariables(), meshNames);

  // remove 1D meshes that were already output by writePolyDataFile
  std::set<std::string> meshesToOutput;
  std::set_difference(meshNames.begin(), meshNames.end(), combined1DMeshes.begin(), combined1DMeshes.end(),
                      std::inserter(meshesToOutput, meshesToOutput.end()));

  // loop over meshes and create a paraview file for each
  for (std::string meshName : meshesToOutput)
  {
    // setup name of file
    std::stringstream filenameStart;
    if (meshesToOutput.size() == 1)
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
    fieldVariable.setRepresentationGlobal();
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
      stringData = Paraview::encodeBase64Float(values.begin(), values.end());
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
        << "type=\"Float32\" "
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
        << "type=\"Float32\" "
        << "NumberOfComponents=\"1\" ";

    std::string stringData;

    // get own rank no
    int ownRankNoCommWorld = 0;
    MPIUtility::handleReturnValue(MPI_Comm_rank(MPI_COMM_WORLD, &ownRankNoCommWorld));

    const node_no_t nNodesLocal = geometryField.functionSpace()->meshPartition()->nNodesLocalWithGhosts();

    std::vector<double> values(nNodesLocal, (double)ownRankNoCommWorld);

    if (binaryOutput)
    {
      stringData = Paraview::encodeBase64Float(values.begin(), values.end());
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

template <typename Iter>
std::string Paraview::encodeBase64Float(Iter iterBegin, Iter iterEnd, bool withEncodedSizePrefix)
{
  // encode as Paraview Float32
  assert(sizeof(float) == 4);

  int dataLength = std::distance(iterBegin, iterEnd)*sizeof(float);
  int rawLength = dataLength;
  int dataStartPos = 0;

  if (withEncodedSizePrefix)
  {
    rawLength += 4;
    dataStartPos = 4;
  }

  int encodedLength = Base64::EncodedLength(rawLength);

  std::vector<char> raw(rawLength, char(0));
  // loop over vector entries and add bytes to raw buffer
  int i = 0;
  for (Iter iter = iterBegin; iter != iterEnd; iter++, i++)
  {
    union {
      float d;    // note, this is a float, not a double, to get 4 bytes
      char c[4];
    };
    d = (float)(*iter);
    memcpy(raw.data() + dataStartPos + i*sizeof(float), c, 4);
  }

  // prepend number of bytes as uint32
  if (withEncodedSizePrefix)
  {
    union {
      uint32_t i;
      char c[4];
    };
    i = rawLength;

    memcpy(raw.data(), c, 4);
  }

  std::vector<char> encoded(encodedLength+1, '\0');

  bool success = Base64::Encode(raw.data(), rawLength, encoded.data(), encodedLength);
  if (!success)
    LOG(WARNING) << "Base64 encoding failed";

  encoded[encodedLength] = '\0';

  return std::string(encoded.begin(), encoded.end());
}

template <typename Iter>
std::string Paraview::encodeBase64Int(Iter iterBegin, Iter iterEnd, bool withEncodedSizePrefix)
{
  // encode as Paraview Int32
  assert(sizeof(int) == 4);

  int dataLength = std::distance(iterBegin, iterEnd)*sizeof(int);
  int rawLength = dataLength;
  int dataStartPos = 0;

  if (withEncodedSizePrefix)
  {
    rawLength += 4;
    dataStartPos = 4;
  }

  int encodedLength = Base64::EncodedLength(rawLength);

  std::vector<char> raw(rawLength, char(0));
  // loop over vector entries and add bytes to raw buffer
  int i = 0;
  for (Iter iter = iterBegin; iter != iterEnd; iter++, i++)
  {
    union {
      int integer;
      char c[4];
    };
    integer = (int)(*iter);
    memcpy(raw.data() + dataStartPos + i*sizeof(float), c, 4);
  }

  // prepend number of bytes as uint32
  if (withEncodedSizePrefix)
  {
    union {
      uint32_t i;
      char c[4];
    };
    i = dataLength;

    memcpy(raw.data(), c, 4);
  }

  std::vector<char> encoded(encodedLength+1, '\0');

  bool success = Base64::Encode(raw.data(), rawLength, (char *)encoded.data(), encodedLength);
  if (!success)
    LOG(WARNING) << "Base64 encoding failed";


  encoded[encodedLength] = '\0';

  return std::string(encoded.begin(), encoded.end());
}


} // namespace
