#include "output_writer/paraview/paraview.h"

#include <fstream>
#include <iostream>
#include <thread>
#include <chrono>
#include <cstdio>  // remove

#include "easylogging++.h"
#include "base64.h"

#include "output_writer/paraview/loop_output.h"
#include "output_writer/paraview/loop_collect_mesh_properties.h"
#include "output_writer/paraview/loop_get_nodal_values.h"
#include "output_writer/paraview/loop_get_geometry_field_nodal_values.h"
#include "output_writer/paraview/poly_data_properties_for_mesh.h"
#include "control/diagnostic_tool/performance_measurement.h"

namespace OutputWriter
{

template<typename FieldVariablesForOutputWriterType>
void Paraview::writePolyDataFile(const FieldVariablesForOutputWriterType &fieldVariables, std::set<std::string> &meshNames)
{
  // output a *.vtp file which contains 1D meshes, if there are any

  bool meshPropertiesInitialized = !meshPropertiesPolyDataFile_.empty();
  std::vector<std::string> meshNamesVector;

  if (!meshPropertiesInitialized)
  {
    Control::PerformanceMeasurement::start("durationParaview1DInit");

    // collect the size data that is needed to compute offsets for parallel file output
    ParaviewLoopOverTuple::loopCollectMeshProperties<FieldVariablesForOutputWriterType>(fieldVariables, meshPropertiesPolyDataFile_, meshNamesVector);

    Control::PerformanceMeasurement::stop("durationParaview1DInit");
  }

  VLOG(1) << "writePolyDataFile on rankSubset_: " << *this->rankSubset_;
  assert(this->rankSubset_);

  VLOG(1) << "meshPropertiesPolyDataFile_: " << meshPropertiesPolyDataFile_;
  /*
  struct PolyDataPropertiesForMesh
  {
    int dimensionality;    //< D=1: object is a VTK "Line", D=2, D=3: object is a VTK "Poly"
    global_no_t nPointsLocal;   //< the number of points needed for representing the mesh, local value of rank
    global_no_t nCellsLocal;    //< the number of VTK "cells", i.e. "Lines" or "Polys", which is the opendihu number of "elements", local value of rank
    global_no_t nPointsGlobal;   //< the number of points needed for representing the mesh, global value of all rank
    global_no_t nCellsGlobal;    //< the number of VTK "cells", i.e. "Lines" or "Polys", which is the opendihu number of "elements", global value of all ranks

    std::vector<PolyDataPropertiesForMesh::DataArrayName> pointDataArrays;   //< <name,nComponents> of PointData DataArray elements
  };
*/

  /* one VTKPiece is the XML element that will be output as <Piece></Piece>. It is created from one or multiple opendihu meshes
   */
  /*
  struct VTKPiece
  {
    std::set<std::string> meshNamesCombinedMeshes;   //< the meshNames of the combined meshes, or only one meshName if it is not a merged mesh
    PolyDataPropertiesForMesh properties;   //< the properties of the merged mesh

    std::string firstScalarName;   //< name of the first scalar field variable of the mesh
    std::string firstVectorName;   //< name of the first non-scalar field variable of the mesh

    //! constructor, initialize nPoints and nCells to 0
    VTKPiece()
    {
      properties.nPointsLocal = 0;
      properties.nCellsLocal = 0;
      properties.nPointsGlobal = 0;
      properties.nCellsGlobal = 0;
      properties.dimensionality = 0;
    }

    //! assign the correct values to firstScalarName and firstVectorName, only if properties has been set
    void setVTKValues()
    {
      // set values for firstScalarName and firstVectorName from the values in pointDataArrays
      for (auto pointDataArray : properties.pointDataArrays)
      {
        if (firstScalarName == "" && pointDataArray.second == 1)
          firstScalarName = pointDataArray.first;
        if (firstVectorName == "" && pointDataArray.second != 1)
          firstVectorName = pointDataArray.first;
      }
    }
  } vtkPiece;
  */

  if (!meshPropertiesInitialized)
  {
    Control::PerformanceMeasurement::start("durationParaview1DInit");

    // parse the collected properties of the meshes that will be output to the file
    for (std::map<std::string,PolyDataPropertiesForMesh>::iterator iter = meshPropertiesPolyDataFile_.begin(); iter != meshPropertiesPolyDataFile_.end(); iter++)
    {
      std::string meshName = iter->first;

      // do not combine meshes other than 1D meshes
      if (iter->second.dimensionality != 1)
        continue;

      // check if this mesh should be combined with other meshes
      bool combineMesh = true;

      // check if mesh can be merged into previous meshes
      if (!vtkPiece1D_.properties.pointDataArrays.empty())   // if properties are already assigned by an earlier mesh
      {
        if (vtkPiece1D_.properties.pointDataArrays.size() != iter->second.pointDataArrays.size())
        {
          LOG(DEBUG) << "Mesh " << meshName << " cannot be combined with " << vtkPiece1D_.meshNamesCombinedMeshes << ". Number of field variables mismatches for "
            << meshName << " (is " << iter->second.pointDataArrays.size() << " instead of " << vtkPiece1D_.properties.pointDataArrays.size() << ")";
          combineMesh = false;
        }
        else
        {
          for (int j = 0; j < iter->second.pointDataArrays.size(); j++)
          {
            if (vtkPiece1D_.properties.pointDataArrays[j].name != iter->second.pointDataArrays[j].name)  // if the name of the jth field variable is different
            {
              LOG(DEBUG) << "Mesh " << meshName << " cannot be combined with " << vtkPiece1D_.meshNamesCombinedMeshes << ". Field variable names mismatch for "
                << meshName << " (there is \"" << vtkPiece1D_.properties.pointDataArrays[j].name << "\" instead of \"" << iter->second.pointDataArrays[j].name << "\")";
              combineMesh = false;
            }
          }
        }

        if (combineMesh)
        {
          VLOG(1) << "Combine mesh " << meshName << " with " << vtkPiece1D_.meshNamesCombinedMeshes
            << ", add " << iter->second.nPointsLocal << " points, " << iter->second.nCellsLocal << " elements to "
            << vtkPiece1D_.properties.nPointsLocal << " points, " << vtkPiece1D_.properties.nCellsLocal << " elements";

          vtkPiece1D_.properties.nPointsLocal += iter->second.nPointsLocal;
          vtkPiece1D_.properties.nCellsLocal += iter->second.nCellsLocal;

          vtkPiece1D_.properties.nPointsGlobal += iter->second.nPointsGlobal;
          vtkPiece1D_.properties.nCellsGlobal += iter->second.nCellsGlobal;
          vtkPiece1D_.setVTKValues();
        }
      }
      else
      {
        VLOG(1) << "this is the first 1D mesh";

        // properties are not yet assigned
        vtkPiece1D_.properties = iter->second; // store properties
        vtkPiece1D_.setVTKValues();
      }

      VLOG(1) << "combineMesh: " << combineMesh;
      if (combineMesh)
      {
        vtkPiece1D_.meshNamesCombinedMeshes.insert(meshName);
      }
    }

    LOG(DEBUG) << "vtkPiece1D_: meshNamesCombinedMeshes: " << vtkPiece1D_.meshNamesCombinedMeshes << ", properties: " << vtkPiece1D_.properties
      << ", firstScalarName: " << vtkPiece1D_.firstScalarName << ", firstVectorName: " << vtkPiece1D_.firstVectorName;

    Control::PerformanceMeasurement::stop("durationParaview1DInit");
  }

  meshNames = vtkPiece1D_.meshNamesCombinedMeshes;

  // if there are no 1D meshes, return
  if (meshNames.empty())
    return;

  if (!meshPropertiesInitialized)
  {
    // add field variable "partitioning" with 1 component
    PolyDataPropertiesForMesh::DataArrayName dataArrayName;
    dataArrayName.name = "partitioning";
    dataArrayName.nComponents = 1;
    dataArrayName.componentNames = std::vector<std::string>(1,"rankNo");

    vtkPiece1D_.properties.pointDataArrays.push_back(dataArrayName);
  }

  // determine filename, broadcast from rank 0
  std::stringstream filename;
  filename << this->filenameBaseWithNo_ << ".vtp";
  int filenameLength = filename.str().length();

  // broadcast length of filename
  MPIUtility::handleReturnValue(MPI_Bcast(&filenameLength, 1, MPI_INT, 0, this->rankSubset_->mpiCommunicator()), "MPI_Bcast (3)");

  std::vector<char> receiveBuffer(filenameLength+1, char(0));
  strcpy(receiveBuffer.data(), filename.str().c_str());
  MPIUtility::handleReturnValue(MPI_Bcast(receiveBuffer.data(), filenameLength, MPI_CHAR, 0, this->rankSubset_->mpiCommunicator()), "MPI_Bcast (4)");

  std::string filenameStr(receiveBuffer.begin(), receiveBuffer.end());

  // remove file if it exists, synchronization afterwards by MPI calls, that is why the remove call is already here
  assert(this->rankSubset_);
  int ownRankNo = this->rankSubset_->ownRankNo();
  if (ownRankNo == 0)
  {
    // open file to ensure that directory exists and file is writable
    std::ofstream file;
    Generic::openFile(file, filenameStr);

    // close and delete file
    file.close();
    std::remove(filenameStr.c_str());
  }

  if (!meshPropertiesInitialized)
  {
    // exchange information about offset in terms of nCells and nPoints
    nCellsPreviousRanks1D_ = 0;
    nPointsPreviousRanks1D_ = 0;
    nPointsGlobal1D_ = 0;
    nLinesGlobal1D_ = 0;

    Control::PerformanceMeasurement::start("durationParaview1DInit");
    Control::PerformanceMeasurement::start("durationParaview1DReduction");
    MPIUtility::handleReturnValue(MPI_Exscan(&vtkPiece1D_.properties.nCellsLocal, &nCellsPreviousRanks1D_, 1, MPI_INT, MPI_SUM, this->rankSubset_->mpiCommunicator()), "MPI_Exscan");
    MPIUtility::handleReturnValue(MPI_Exscan(&vtkPiece1D_.properties.nPointsLocal, &nPointsPreviousRanks1D_, 1, MPI_INT, MPI_SUM, this->rankSubset_->mpiCommunicator()), "MPI_Exscan");
    MPIUtility::handleReturnValue(MPI_Reduce(&vtkPiece1D_.properties.nPointsLocal, &nPointsGlobal1D_, 1, MPI_INT, MPI_SUM, 0, this->rankSubset_->mpiCommunicator()), "MPI_Reduce");
    MPIUtility::handleReturnValue(MPI_Reduce(&vtkPiece1D_.properties.nCellsLocal, &nLinesGlobal1D_, 1, MPI_INT, MPI_SUM, 0, this->rankSubset_->mpiCommunicator()), "MPI_Reduce");
    Control::PerformanceMeasurement::stop("durationParaview1DReduction");
    Control::PerformanceMeasurement::stop("durationParaview1DInit");
  }

  // get local data values
  // setup connectivity array
  std::vector<int> connectivityValues(2*vtkPiece1D_.properties.nCellsLocal);
  for (int i = 0; i < vtkPiece1D_.properties.nCellsLocal; i++)
  {
    connectivityValues[2*i + 0] = nPointsPreviousRanks1D_ + i;
    connectivityValues[2*i + 1] = nPointsPreviousRanks1D_ + i+1;
  }

  // setup offset array
  std::vector<int> offsetValues(vtkPiece1D_.properties.nCellsLocal);
  for (int i = 0; i < vtkPiece1D_.properties.nCellsLocal; i++)
  {
    offsetValues[i] = 2*nCellsPreviousRanks1D_ + 2*i + 1;
  }

  // collect all data for the field variables, organized by field variable names
  std::map<std::string, std::vector<double>> fieldVariableValues;
  ParaviewLoopOverTuple::loopGetNodalValues<FieldVariablesForOutputWriterType>(fieldVariables, vtkPiece1D_.meshNamesCombinedMeshes, fieldVariableValues);

  assert (!fieldVariableValues.empty());
  fieldVariableValues["partitioning"].resize(vtkPiece1D_.properties.nPointsLocal, (double)this->rankSubset_->ownRankNo());

  // if next assertion will fail, output why for debugging
  if (fieldVariableValues.size() != vtkPiece1D_.properties.pointDataArrays.size())
  {
    LOG(DEBUG) << "n field variable values: " << fieldVariableValues.size() << ", n point data arrays: "
      << vtkPiece1D_.properties.pointDataArrays.size();
    LOG(DEBUG) << "vtkPiece1D_.meshNamesCombinedMeshes: " << vtkPiece1D_.meshNamesCombinedMeshes;
    std::stringstream pointDataArraysNames;
    for (int i = 0; i < vtkPiece1D_.properties.pointDataArrays.size(); i++)
    {
      pointDataArraysNames << vtkPiece1D_.properties.pointDataArrays[i].name << " ";
    }
    LOG(DEBUG) << "pointDataArraysNames: " <<  pointDataArraysNames.str();
  }

  assert(fieldVariableValues.size() == vtkPiece1D_.properties.pointDataArrays.size());

#ifndef NDEBUG
  LOG(DEBUG) << "fieldVariableValues: ";
  for (std::map<std::string, std::vector<double>>::iterator iter = fieldVariableValues.begin(); iter != fieldVariableValues.end(); iter++)
  {
    LOG(DEBUG) << iter->first;
  }
#endif

  // check if field variable names have changed since last initialization
  for (std::vector<PolyDataPropertiesForMesh::DataArrayName>::iterator pointDataArrayIter = vtkPiece1D_.properties.pointDataArrays.begin();
       pointDataArrayIter != vtkPiece1D_.properties.pointDataArrays.end(); pointDataArrayIter++)
  {
    LOG(DEBUG) << "  field variable \"" << pointDataArrayIter->name << "\".";

    // if there is a field variable with a name that was not present when vtkPiece1D_ was created
    if (fieldVariableValues.find(pointDataArrayIter->name) == fieldVariableValues.end())
    {
      LOG(DEBUG) << "Field variable names have changed, reinitialize Paraview output writer.";

      // reset now old variables
      meshPropertiesInitialized = false;
      meshPropertiesPolyDataFile_.clear();
      vtkPiece1D_ = VTKPiece();

      // recursively call this method
      writePolyDataFile(fieldVariables, meshNames);

      return;
    }
  }

  // collect all data for the geometry field variable
  std::vector<double> geometryFieldValues;
  ParaviewLoopOverTuple::loopGetGeometryFieldNodalValues<FieldVariablesForOutputWriterType>(fieldVariables, vtkPiece1D_.meshNamesCombinedMeshes, geometryFieldValues);

  // only continue if there is data to reduce
  if (vtkPiece1D_.meshNamesCombinedMeshes.empty())
  {
    LOG(ERROR) << "There are no 1D meshes that could be combined, but Paraview output with combineFiles=True was specified. \n(This only works for 1D meshes.)";
  }

  LOG(DEBUG) << "Combined mesh from " << vtkPiece1D_.meshNamesCombinedMeshes;

  int nOutputFileParts = 4 + vtkPiece1D_.properties.pointDataArrays.size();

  // transform current time to string
  std::vector<double> time(1, this->currentTime_);
  std::string stringTime;
  if (binaryOutput_)
  {
    stringTime = Paraview::encodeBase64Float(time.begin(), time.end());
  }
  else
  {
    stringTime = Paraview::convertToAscii(time, fixedFormat_);
  }

  // create the basic structure of the output file
  std::vector<std::stringstream> outputFileParts(nOutputFileParts);
  int outputFilePartNo = 0;
  outputFileParts[outputFilePartNo] << "<?xml version=\"1.0\"?>" << std::endl
    << "<!-- " << DihuContext::versionText() << " " << DihuContext::metaText()
    << ", currentTime: " << this->currentTime_ << ", timeStepNo: " << this->timeStepNo_ << " -->" << std::endl
    << "<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\">" << std::endl    // intel cpus are LittleEndian
    << std::string(1, '\t') << "<PolyData>" << std::endl
    << std::string(2, '\t') << "<FieldData>" << std::endl
    << std::string(3, '\t') << "<DataArray type=\"Float32\" Name=\"Time\" NumberOfTuples=\"1\" format=\"" << (binaryOutput_? "binary" : "ascii")
    << "\" >" << std::endl
    << std::string(4, '\t') << stringTime << std::endl
    << std::string(3, '\t') << "</DataArray>" << std::endl
    << std::string(2, '\t') << "</FieldData>" << std::endl;

  outputFileParts[outputFilePartNo] << std::string(2, '\t') << "<Piece NumberOfPoints=\"" << nPointsGlobal1D_ << "\" NumberOfVerts=\"0\" "
    << "NumberOfLines=\"" << nLinesGlobal1D_ << "\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">" << std::endl
    << std::string(3, '\t') << "<PointData";

  if (vtkPiece1D_.firstScalarName != "")
  {
    outputFileParts[outputFilePartNo] << " Scalars=\"" << vtkPiece1D_.firstScalarName << "\"";
  }
  if (vtkPiece1D_.firstVectorName != "")
  {
    outputFileParts[outputFilePartNo] << " Vectors=\"" << vtkPiece1D_.firstVectorName << "\"";
  }
  outputFileParts[outputFilePartNo] << ">" << std::endl;

  // loop over field variables (PointData)
  for (std::vector<PolyDataPropertiesForMesh::DataArrayName>::iterator pointDataArrayIter = vtkPiece1D_.properties.pointDataArrays.begin(); pointDataArrayIter != vtkPiece1D_.properties.pointDataArrays.end(); pointDataArrayIter++)
  {
    // set up string for component names
    std::stringstream componentNames;
    bool isComponentNamesSet = false;
    for (int componentNo = 0; componentNo < pointDataArrayIter->nComponents; componentNo++)
    {
      // get component name of the current component
      std::string componentName = pointDataArrayIter->componentNames[componentNo];

      componentNames << "ComponentName" << componentNo << "=\"" << componentName << "\" ";

      // check if it is equal to simply the componentNo, if no component name was explicitly specify they have names "0", "1", ...
      std::stringstream trivialComponentName;
      trivialComponentName << componentNo;

      if (componentName != trivialComponentName.str())
        isComponentNamesSet = true;
    }

    // if there were no real component names set, do not include them in the VTK file
    if (!isComponentNamesSet)
      componentNames.str("");

    // write normal data element
    outputFileParts[outputFilePartNo] << std::string(4, '\t') << "<DataArray "
        << "Name=\"" << pointDataArrayIter->name << "\" "
        << "type=\"" << (pointDataArrayIter->name == "partitioning"? "Int32" : "Float32") << "\" "
        << "NumberOfComponents=\"" << pointDataArrayIter->nComponents << "\" "
        << componentNames.str()
        << "format=\"" << (binaryOutput_? "binary" : "ascii")
        << "\" >" << std::endl << std::string(5, '\t');

    // at this point the data of the field variable is missing
    outputFilePartNo++;

    outputFileParts[outputFilePartNo] << std::endl << std::string(4, '\t') << "</DataArray>" << std::endl;
  }

  outputFileParts[outputFilePartNo] << std::string(3, '\t') << "</PointData>" << std::endl
    << std::string(3, '\t') << "<CellData>" << std::endl
    << std::string(3, '\t') << "</CellData>" << std::endl
    << std::string(3, '\t') << "<Points>" << std::endl
    << std::string(4, '\t') << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"" << (binaryOutput_? "binary" : "ascii")
    << "\" >" << std::endl << std::string(5, '\t');

  // at this point the data of points (geometry field) is missing
  outputFilePartNo++;

  outputFileParts[outputFilePartNo] << std::endl << std::string(4, '\t') << "</DataArray>" << std::endl
    << std::string(3, '\t') << "</Points>" << std::endl
    << std::string(3, '\t') << "<Verts></Verts>" << std::endl
    << std::string(3, '\t') << "<Lines>" << std::endl
    << std::string(4, '\t') << "<DataArray Name=\"connectivity\" type=\"Int32\" "
    << (binaryOutput_? "format=\"binary\"" : "format=\"ascii\"") << ">" << std::endl << std::string(5, '\t');

  // at this point the the structural information of the lines (connectivity) is missing
  outputFilePartNo++;

  outputFileParts[outputFilePartNo]
    << std::endl << std::string(4, '\t') << "</DataArray>" << std::endl
    << std::string(4, '\t') << "<DataArray Name=\"offsets\" type=\"Int32\" "
    << (binaryOutput_? "format=\"binary\"" : "format=\"ascii\"") << ">" << std::endl << std::string(5, '\t');

  // at this point the offset array will be written to the file
  outputFilePartNo++;

  outputFileParts[outputFilePartNo]
    << std::endl << std::string(4, '\t') << "</DataArray>" << std::endl
    << std::string(3, '\t') << "</Lines>" << std::endl
    << std::string(3, '\t') << "<Strips></Strips>" << std::endl
    << std::string(3, '\t') << "<Polys></Polys>" << std::endl
    << std::string(2, '\t') << "</Piece>" << std::endl
    << std::string(1, '\t') << "</PolyData>" << std::endl
    << "</VTKFile>" << std::endl;

  assert(outputFilePartNo+1 == nOutputFileParts);

  // loop over output file parts and collect the missing data for the own rank

  VLOG(1) << "outputFileParts:";

  for (std::vector<std::stringstream>::iterator iter = outputFileParts.begin(); iter != outputFileParts.end(); iter++)
  {
    VLOG(1) << "  " << iter->str();
  }

  LOG(DEBUG) << "open MPI file \"" << filenameStr << "\".";

  // open file
  MPI_File fileHandle;
  MPIUtility::handleReturnValue(MPI_File_open(this->rankSubset_->mpiCommunicator(), filenameStr.c_str(),
                                              //MPI_MODE_WRONLY | MPI_MODE_CREATE | MPI_MODE_UNIQUE_OPEN,
                                              MPI_MODE_WRONLY | MPI_MODE_CREATE,
                                              MPI_INFO_NULL, &fileHandle), "MPI_File_open");

  Control::PerformanceMeasurement::start("durationParaview1DWrite");

  // write beginning of file on rank 0
  outputFilePartNo = 0;

  writeAsciiDataShared(fileHandle, ownRankNo, outputFileParts[outputFilePartNo].str());
  outputFilePartNo++;

  VLOG(1) << "get current shared file position";

  // get current file position
  MPI_Offset currentFilePosition = 0;
  MPIUtility::handleReturnValue(MPI_File_get_position_shared(fileHandle, &currentFilePosition), "MPI_File_get_position_shared");
  LOG(DEBUG) << "current shared file position: " << currentFilePosition;

  // write field variables
  // loop over field variables
  int fieldVariableNo = 0;
  for (std::vector<PolyDataPropertiesForMesh::DataArrayName>::iterator pointDataArrayIter = vtkPiece1D_.properties.pointDataArrays.begin();
       pointDataArrayIter != vtkPiece1D_.properties.pointDataArrays.end(); pointDataArrayIter++, fieldVariableNo++)
  {
    assert(fieldVariableValues.find(pointDataArrayIter->name) != fieldVariableValues.end());

    // write values
    bool writeFloatsAsInt = pointDataArrayIter->name == "partitioning";    // for partitioning, convert float values to integer values for output
    writeCombinedValuesVector(fileHandle, ownRankNo, fieldVariableValues[pointDataArrayIter->name], fieldVariableNo, writeFloatsAsInt);

    // write next xml constructs
    writeAsciiDataShared(fileHandle, ownRankNo, outputFileParts[outputFilePartNo].str());
    outputFilePartNo++;
  }

  // write geometry field data
  writeCombinedValuesVector(fileHandle, ownRankNo, geometryFieldValues, fieldVariableNo++);

  // write next xml constructs
  writeAsciiDataShared(fileHandle, ownRankNo, outputFileParts[outputFilePartNo].str());
  outputFilePartNo++;

  // write connectivity values
  writeCombinedValuesVector(fileHandle, ownRankNo, connectivityValues, fieldVariableNo++);

  // write next xml constructs
  writeAsciiDataShared(fileHandle, ownRankNo, outputFileParts[outputFilePartNo].str());
  outputFilePartNo++;

  // write offset values
  writeCombinedValuesVector(fileHandle, ownRankNo, offsetValues, fieldVariableNo++);

  // write next xml constructs
  writeAsciiDataShared(fileHandle, ownRankNo, outputFileParts[outputFilePartNo].str());

  /*
    int array_of_sizes[1];
    array_of_sizes[0]=numProcs;
    int array_of_subsizes[1];
    array_of_subsizes[0]=1;
    int array_of_starts[1];
    array_of_starts[0]=myId;


    MPI_Datatype accessPattern;
    MPI_Type_create_subarray(1,array_of_sizes, array_of_subsizes, array_of_starts, MPI_ORDER_C, MPI_BYTE, &accessPattern);
    MPI_Type_commit(&accessPattern);

    MPI_File_set_view(fh, 0, MPI_BYTE, accessPattern, "native", MPI_INFO_NULL);
    MPI_File_write(fh, v, size, MPI_BYTE, MPI_STATUS_IGNORE);
  */

  Control::PerformanceMeasurement::stop("durationParaview1DWrite");

  MPIUtility::handleReturnValue(MPI_File_close(&fileHandle), "MPI_File_close");

  // register file at SeriesWriter to be included in the "*.vtk.series" JSON file
  if (ownRankNo == 0)
  {
    Paraview::seriesWriter().registerNewFile(filenameStr, this->currentTime_);
  }
}

} // namespace
