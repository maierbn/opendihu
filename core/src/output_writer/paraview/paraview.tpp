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
#include "control/performance_measurement.h"

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

  Control::PerformanceMeasurement::start("durationParaviewOutput");

  std::set<std::string> combined1DMeshes;
  std::set<std::string> combined2DMeshes;
  std::set<std::string> combined3DMeshes;

  if (combineFiles_)
  {
    Control::PerformanceMeasurement::start("durationParaview1D");

    // create a PolyData file that combines all 1D meshes into one file
    writePolyDataFile<typename DataType::FieldVariablesForOutputWriter>(data.getFieldVariablesForOutputWriter(), combined1DMeshes);

    Control::PerformanceMeasurement::stop("durationParaview1D");
    Control::PerformanceMeasurement::start("durationParaview3D");

    // create an UnstructuredMesh file that combines all 3D meshes into one file
    writeCombinedUnstructuredGridFile<typename DataType::FieldVariablesForOutputWriter>(data.getFieldVariablesForOutputWriter(), combined3DMeshes, true);

    Control::PerformanceMeasurement::stop("durationParaview3D");
    Control::PerformanceMeasurement::start("durationParaview2D");

    // create an UnstructuredMesh file that combines all 2D meshes into one file
    writeCombinedUnstructuredGridFile<typename DataType::FieldVariablesForOutputWriter>(data.getFieldVariablesForOutputWriter(), combined2DMeshes, false);

    Control::PerformanceMeasurement::stop("durationParaview2D");
  }

  // output normal files, parallel or if combineFiles_, only the 2D and 3D meshes, combined

  // collect all available meshes
  std::set<std::string> meshNames;
  LoopOverTuple::loopCollectMeshNames<typename DataType::FieldVariablesForOutputWriter>(data.getFieldVariablesForOutputWriter(), meshNames);

  // remove 1D meshes that were already output by writePolyDataFile
  std::set<std::string> meshesWithout1D;
  std::set_difference(meshNames.begin(), meshNames.end(), combined1DMeshes.begin(), combined1DMeshes.end(),
                      std::inserter(meshesWithout1D, meshesWithout1D.end()));

  // remove 3D meshes that were already output by writeCombinedUnstructuredGridFile
  std::set<std::string> meshesWithout1D3D;
  std::set_difference(meshesWithout1D.begin(), meshesWithout1D.end(), combined3DMeshes.begin(), combined3DMeshes.end(),
                      std::inserter(meshesWithout1D3D, meshesWithout1D3D.end()));

  // remove 2D meshes that were already output by writeCombinedUnstructuredGridFile
  std::set<std::string> meshesToOutput;
  std::set_difference(meshesWithout1D3D.begin(), meshesWithout1D3D.end(), combined2DMeshes.begin(), combined2DMeshes.end(),
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
    ParaviewLoopOverTuple::loopOutput(data.getFieldVariablesForOutputWriter(), data.getFieldVariablesForOutputWriter(), meshName, filenameStart.str(), specificSettings_);
  }

  Control::PerformanceMeasurement::stop("durationParaviewOutput");
}

template<typename FieldVariableType>
void Paraview::writeParaviewFieldVariable(FieldVariableType &fieldVariable,
                                          std::ofstream &file, bool binaryOutput, bool fixedFormat, bool onlyParallelDatasetElement)
{
  // here we have the type of the mesh with meshName (which is typedef to FunctionSpace)
  //typedef typename FieldVariableType::FunctionSpace FunctionSpace;

  // paraview does not correctly handle 2-component output data, so set number to 3
  int nComponentsParaview = fieldVariable.nComponents();
  if (nComponentsParaview == 2)
    nComponentsParaview = 3;

  // if only the "parallel dataset element" stub which is needed in the master files, should be written
  if (onlyParallelDatasetElement)
  {

    file << std::string(3, '\t') << "<PDataArray "
        << "Name=\"" << fieldVariable.name() << "\" "
        << "type=\"Float32\" "
        << "NumberOfComponents=\"" << nComponentsParaview << "\" ";

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
        << "NumberOfComponents=\"" << nComponentsParaview << "\" ";

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
    values.reserve(componentValues[0].size()*nComponentsParaview);

    // copy values in consecutive order (x y z x y z) to output
    for (int i = 0; i < componentValues[0].size(); i++)
    {
      for (int componentNo = 0; componentNo < nComponentsParaview; componentNo++)
      {
        if (nComponents == 2 && nComponentsParaview == 3 && componentNo == 2)
        {
          values.push_back(0.0);
        }
        else
        {
          values.push_back(componentValues[componentNo][i]);
        }
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

    const node_no_t nNodesLocal = geometryField.functionSpace()->meshPartition()->nNodesLocalWithGhosts();

    std::vector<int32_t> values(nNodesLocal, (int32_t)DihuContext::ownRankNoCommWorld());

    if (binaryOutput)
    {
      stringData = Paraview::encodeBase64Int32(values.begin(), values.end());
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

  std::vector<char> encoded(encodedLength, '\0');

  bool success = Base64::Encode(raw.data(), rawLength, encoded.data(), encodedLength);
  if (!success)
    LOG(WARNING) << "Base64 encoding failed";

  return std::string(encoded.begin(), encoded.end());
}

template <typename Iter>
std::string Paraview::encodeBase64Int32(Iter iterBegin, Iter iterEnd, bool withEncodedSizePrefix)
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

  std::vector<char> encoded(encodedLength, '\0');

  bool success = Base64::Encode(raw.data(), rawLength, (char *)encoded.data(), encodedLength);
  if (!success)
    LOG(WARNING) << "Base64 encoding failed";

  return std::string(encoded.begin(), encoded.end());
}

template <typename Iter>
std::string Paraview::encodeBase64UInt8(Iter iterBegin, Iter iterEnd, bool withEncodedSizePrefix)
{
  // encode as Paraview UInt8

  int dataLength = std::distance(iterBegin, iterEnd);
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
    raw[dataStartPos + i] = char(*iter);
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

  std::vector<char> encoded(encodedLength, '\0');

  bool success = Base64::Encode(raw.data(), rawLength, (char *)encoded.data(), encodedLength);
  if (!success)
    LOG(WARNING) << "Base64 encoding failed";

  return std::string(encoded.begin(), encoded.end());
}


} // namespace
