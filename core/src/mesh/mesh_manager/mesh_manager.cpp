#include "mesh/mesh_manager/mesh_manager.h"

#include "function_space/function_space.h"
#include "mesh/structured_regular_fixed.h"
#include "mesh/unstructured_deformable.h"

namespace Mesh
{

Manager::Manager(PythonConfig specificSettings) :
  partitionManager_(nullptr), specificSettings_(specificSettings), numberAnonymousMeshes_(0)
{
  LOG(TRACE) << "MeshManager constructor";
  storePreconfiguredMeshes();
}

void Manager::setPartitionManager(std::shared_ptr<Partition::Manager> partitionManager)
{
  partitionManager_ = partitionManager;
}

void Manager::storePreconfiguredMeshes()
{
  LOG(TRACE) << "MeshManager::storePreconfiguredMeshes";
  if (this->specificSettings_.pyObject())
  {
    std::string keyString("Meshes");
    if (specificSettings_.hasKey("Meshes"))
    {

      // loop over entries of python dict "Meshes"
      std::pair<std::string, PyObject *> dictItem
        = specificSettings_.getOptionDictBegin<std::string, PyObject *>(keyString);

      for (; !specificSettings_.getOptionDictEnd(keyString);
          specificSettings_.getOptionDictNext<std::string, PyObject *>(keyString, dictItem))
      {
        // the key is the mesh name
        std::string key = dictItem.first;

        // the value is the config for the mesh
        PyObject *value = dictItem.second;

        if (value == NULL)
        {
          LOG(WARNING) << "Could not extract dict for Mesh \"" << key << "\".";
        }
        else if (!PyDict_Check(value))
        {
          LOG(WARNING) << "Value for mesh with name \"" << key << "\" should be a dict.";
        }
        else
        {
          VLOG(1) << "Store mesh configuration with key \"" << key << "\".";
          if (meshConfiguration_.find(key) != meshConfiguration_.end())
          {
            meshConfiguration_.at(key).setPyObject(value);
          }
          else
          {
            meshConfiguration_.insert(std::pair<std::string,PythonConfig>(key, PythonConfig(specificSettings_, "Meshes", key, value)));
          }

          // check if mesh contains node positions as file name and offset
          if (PythonUtility::hasKey(value, "nodePositions"))
          {
            PyObject *nodePositions = PythonUtility::getOptionPyObject(value, "nodePositions", specificSettings_.getStringPath());
            if (nodePositions)
            {
              if (PythonUtility::isTypeList(nodePositions))
              {
                // get first entry of list
                PyObject *item = PythonUtility::getOptionListBegin<PyObject *>(value, "nodePositions", specificSettings_.getStringPath());
                if (!PythonUtility::getOptionListEnd(value, "nodePositions", specificSettings_.getStringPath()))
                {
                  if (PyUnicode_Check(item))
                  {

                    std::string filename = PythonUtility::convertFromPython<std::string>::get(item);
                    PythonUtility::getOptionListNext<PyObject *>(value, "nodePositions", specificSettings_.getStringPath(), item);

                    LOG(DEBUG) << "extracted filename \"" << filename << "\"";

                    if (!PythonUtility::getOptionListEnd(value, "nodePositions", specificSettings_.getStringPath()))
                    {
                      nodePositionsFromFile_[key] = NodePositionsFromFile();
                      nodePositionsFromFile_[key].filename = filename;

                      nodePositionsFromFile_[key].chunks = PythonUtility::convertFromPython<std::vector<std::pair<MPI_Offset,int>>>::get(item);
                    }
                    else
                    {
                      LOG(WARNING) << specificSettings_.getStringPath() << "[\"nodePositions\"]: If nodePositions contains a filename, then a list [[offset,length],[offset,length]...] has to follow.";
                    }
                  }
                }
              }
            }
          }
        }
      }

      loadGeometryFromFile();
    }
    else
    {
      LOG(INFO) << "You have specified the mesh in-line and not under the extra key \"Meshes\". You could do so,"
        " by defining \"Meshes\": {\"<your custom mesh name>\": {<your mesh parameters>}} at the beginning of the"
        " config and \"meshName\": \"<your custom mesh name>\" where you currently have specified the mesh parameters."
        " This is required if you want to use the same mesh for multiple objects.";
    }
  }
}

void Manager::loadGeometryFromFile()
{
  Control::PerformanceMeasurement::start("durationReadGeometry");

  // at this point all ranks in MPI_COMM_WORLD execute collectively
  MPI_Comm mpiCommunicator = MPI_COMM_WORLD;;
  // check if some nodePositions of at least one rank are in a file
  int nodePositionsAreInFileOnRank = -1;   // -1 means there are no nodePositions in a file
  if (!nodePositionsFromFile_.empty())
  {
    nodePositionsAreInFileOnRank = DihuContext::ownRankNoCommWorld();
  }
  int nodePositionsAreInFileOnRankGlobal = 0;
  MPIUtility::handleReturnValue(MPI_Allreduce(&nodePositionsAreInFileOnRank, &nodePositionsAreInFileOnRankGlobal, 1, MPI_INT, MPI_MAX, mpiCommunicator), "MPI_Allreduce");

  LOG(DEBUG) << "nodePositionsAreInFileOnRankGlobal: " << nodePositionsAreInFileOnRankGlobal;

  if (nodePositionsAreInFileOnRankGlobal != -1)
  {
    // now nodePositionsAreInFileOnRankGlobal holds a rank that has a file name
    std::string filename;
    if (DihuContext::ownRankNoCommWorld() == nodePositionsAreInFileOnRankGlobal)
    {
      filename = nodePositionsFromFile_.begin()->second.filename;
      LOG(DEBUG) << "there is a file on this rank, " << filename;
    }
    int filenameLength = filename.length();

    // broadcast length of filename
    MPIUtility::handleReturnValue(MPI_Bcast(&filenameLength, 1, MPI_INT, nodePositionsAreInFileOnRankGlobal, mpiCommunicator), "MPI_Bcast");

    // broadcast filename
    std::vector<char> receiveBuffer(filenameLength+1, char(0));
    strcpy(receiveBuffer.data(), filename.c_str());
    MPIUtility::handleReturnValue(MPI_Bcast(receiveBuffer.data(), filenameLength, MPI_CHAR, nodePositionsAreInFileOnRankGlobal, mpiCommunicator), "MPI_Bcast");

    std::string openFileName(receiveBuffer.begin(), receiveBuffer.end());
    while (openFileName[openFileName.length()-1] == 0)
    {
      openFileName = openFileName.substr(0,openFileName.length()-1);
    }

    LOG(DEBUG) << "openFileName: " << openFileName;

    // MPI-reduce maximum number of times (number of meshes) to read
    int nTimesInFile = 0;
    for (std::map<std::string, NodePositionsFromFile>::iterator iter = nodePositionsFromFile_.begin();
        iter != nodePositionsFromFile_.end(); iter++)
    {
      std::string filename = iter->second.filename;
      LOG(DEBUG) << "check filename [" << filename << "] == [" << openFileName << "] length: " << filename.length() << "," << openFileName.length();
      if (std::string(filename) == std::string(openFileName))
      {
        LOG(DEBUG) << "yes";
        nTimesInFile += iter->second.chunks.size();
      }
    }
    int nTimesInFileGlobal = 0;
    MPIUtility::handleReturnValue(MPI_Allreduce(&nTimesInFile, &nTimesInFileGlobal, 1, MPI_INT, MPI_MAX, mpiCommunicator), "MPI_Allreduce");

    LOG(DEBUG) << "number times that node positions for a mesh are given by a file: local: " << nTimesInFile << ", global: " << nTimesInFileGlobal;

    LOG(INFO) << "Read from file \"" << openFileName << "\", " << nTimesInFileGlobal << " collective chunks.";

    // collectively open the file for reading
    MPI_File fileHandle;
    MPIUtility::handleReturnValue(MPI_File_open(mpiCommunicator, openFileName.c_str(),
                                                //MPI_MODE_WRONLY | MPI_MODE_CREATE | MPI_MODE_UNIQUE_OPEN,
                                                MPI_MODE_RDONLY,
                                                MPI_INFO_NULL, &fileHandle), "MPI_File_open");

    double progress = 0;
    // loop over node positions in file to be read
    int nNodePositionsRead = 0;
    for (std::map<std::string, NodePositionsFromFile>::iterator nodePositionsFromFileIter = nodePositionsFromFile_.begin();
        nodePositionsFromFileIter != nodePositionsFromFile_.end(); nodePositionsFromFileIter++)
    {
      NodePositionsFromFile &nodePositions = nodePositionsFromFileIter->second;
      std::string filename = nodePositionsFromFileIter->second.filename;
      if (filename == openFileName)
      {
        // loop over chunks
        for (std::vector<std::pair<MPI_Offset,int>>::iterator chunksIter = nodePositions.chunks.begin(); chunksIter != nodePositions.chunks.end(); chunksIter++)
        {
          double newProgress = (double)nNodePositionsRead / nTimesInFileGlobal;
          if (int(newProgress*10) != int(progress*10))
          {
            progress = newProgress;
            if (DihuContext::ownRankNoCommWorld() == 0)
            {
              std::cout << "\b\b\b\b" << int(progress*100) << "%" << std::flush;
            }
          }
          progress = newProgress;

          MPI_Offset offset = chunksIter->first;
          int nValues = chunksIter->second*3;

          int oldSize = nodePositions.data.size();
          int newSize = oldSize + nValues;
          nodePositions.data.resize(newSize);

          MPIUtility::handleReturnValue(MPI_File_read_at_all(fileHandle, offset, nodePositions.data.data()+oldSize, nValues, MPI_DOUBLE, MPI_STATUS_IGNORE), "MPI_Read_at_all");

          LOG(DEBUG) << "for mesh \"" << nodePositionsFromFileIter->first << "\" read from file \"" << openFileName << "\" offset " << offset
            << ", nValues (number of values): " << nValues << ", last 3 values: "
            << "[" << nodePositions.data[nodePositions.data.size()-3] << "," << nodePositions.data[nodePositions.data.size()-2] << "," << nodePositions.data[nodePositions.data.size()-1] << "]";

          nNodePositionsRead++;
        }
      }
    }
    LOG(DEBUG) << "collectively read other node positions, only on other ranks (" << nTimesInFileGlobal-nNodePositionsRead << " remaining)";
    for(; nNodePositionsRead < nTimesInFileGlobal; nNodePositionsRead++)
    {
      double newProgress = (double)nNodePositionsRead / nTimesInFileGlobal;
      if (int(newProgress*10) != int(progress*10))
      {
        progress = newProgress;
        if (DihuContext::ownRankNoCommWorld() == 0)
        {
          std::cout << "\b\b\b\b" << int(progress*100) << "%" << std::flush;
        }
      }
      progress = newProgress;

      std::vector<double> readBuffer(1);
      MPIUtility::handleReturnValue(MPI_File_read_at_all(fileHandle, 0, readBuffer.data(), 0, MPI_DOUBLE, MPI_STATUS_IGNORE), "MPI_Read_at_all");
    }
    if (DihuContext::ownRankNoCommWorld() == 0)
    {
      std::cout << "\b\b\b\bdone." << std::endl;
    }

    MPIUtility::handleReturnValue(MPI_File_close(&fileHandle), "MPI_File_close");
  }

  Control::PerformanceMeasurement::stop("durationReadGeometry");
}

std::shared_ptr<FunctionSpace::Generic> Manager::
createGenericFunctionSpace(int nEntries, std::string name)
{
  // constructor is declared in function_space/06_function_space_dofs_nodes.h
  // FunctionSpaceDofsNodes(std::shared_ptr<Partition::Manager> partitionManager, std::array<element_no_t, D> nElements, std::array<double, D> physicalExtent, int inputMeshIsGlobal);

  std::array<element_no_t, 1> nElements({nEntries - 1});
  std::array<double, 1> physicalExtent({0.0});
  std::array<int, 1> nRanksPerCoordinateDirection({1});
  std::shared_ptr<Mesh> mesh = createFunctionSpace<FunctionSpace::Generic>(name, nElements, physicalExtent, nRanksPerCoordinateDirection, false);   // last parameter is that nElements is local number

  return std::static_pointer_cast<FunctionSpace::Generic>(mesh);
}

std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace::Generic,1>> Manager::
createGenericFieldVariable(int nEntries, std::string name)
{
  assert(nEntries > 1);

  // create generic field variable
  std::stringstream meshName;
  meshName << "meshForFieldVariable" << name;
  LOG(DEBUG) << "create generic field variable with " << nEntries << " entries.";
  std::shared_ptr<FunctionSpace::Generic> functionSpace = createGenericFunctionSpace(nEntries, meshName.str());

  // createFieldVariable is declared in function_space/10_function_space_field_variable.h
  //template <int nComponents>
  //std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace<MeshType,BasisFunctionType>,nComponents>> createFieldVariable(std::string name);
  return functionSpace->template createFieldVariable<1>(name);
}

bool Manager::hasFunctionSpace(std::string meshName)
{
  return functionSpaces_.find(meshName) != functionSpaces_.end();
}

//! remove a function space if it exists
void Manager::deleteFunctionSpace(std::string meshName)
{
  if (hasFunctionSpace(meshName))
  {
    LOG(DEBUG) << "deleteFunctionSpace(" << meshName << ")";
    functionSpaces_.erase(meshName);
  }
}

} // namespace
