#include "output_writer/output_surface/output_surface.h"

#include <algorithm>

namespace OutputWriter
{

template<typename Solver>
OutputSurface<Solver>::
OutputSurface(DihuContext context) :
  context_(context["OutputSurface"]), solver_(context_),
  data_(context_), ownRankInvolvedInOutput_(true), timeStepNo_(0), currentTime_(0.0)
{

}

template<typename Solver>
void OutputSurface<Solver>::
advanceTimeSpan()
{
  solver_.advanceTimeSpan();

  LOG(DEBUG) << "OutputSurface: writeOutput, ownRankInvolvedInOutput_: " << ownRankInvolvedInOutput_;
  if (ownRankInvolvedInOutput_)
  {
    outputWriterManager_.writeOutput(data_, timeStepNo_++);
  }

  // write out values at points
  writeSampledPoints();
}

template<typename Solver>
void OutputSurface<Solver>::
initialize()
{
  if (initialized_)
    return;


  // initialize solvers
  solver_.initialize();
  data_.setData(solver_.data());
  data_.initialize();
  ownRankInvolvedInOutput_ = data_.ownRankInvolvedInOutput();

  // initialize output writers
  PythonConfig specificSettings = context_.getPythonConfig();

  // initialize points
  PyObject *samplingPointsPy = specificSettings.getOptionPyObject("samplingPoints");
  sampledPoints_ = PythonUtility::convertFromPython<std::vector<Vec3>>::get(samplingPointsPy);

  if (!sampledPoints_.empty())
    filename_ = specificSettings.getOptionString("filename", "out/sampledPoints.csv");

  initializeSampledPoints();

  LOG(DEBUG) << "OutputSurface: initialize output writers";

  // initialize output writer to use smaller rank subset that only contains the ranks that have parts of the surface
  // if the last argument is not given, by default the common rank subset would be used
  if (ownRankInvolvedInOutput_)
  {
    rankSubset_ = data_.functionSpace()->meshPartition()->rankSubset();
    outputWriterManager_.initialize(context_, specificSettings, rankSubset_);
  }
  initialized_ = true;
}

template<typename Solver>
void OutputSurface<Solver>::
initializeSampledPoints()
{
  if (sampledPoints_.empty())
    return;

  LOG(DEBUG) << "initialize sampled points";

  // get the 2D function spaces
  this->data_.getFunctionSpaces(functionSpaces_);

  // now we have a 2D function space for each face in functionSpaces

  // loop over sampling points and find them in the function spaces
  for (int samplingPointNo = 0; samplingPointNo < sampledPoints_.size(); samplingPointNo++)
  {
    // loop over function spaces for the faces
    for (int functionSpaceNo = 0; functionSpaceNo < functionSpaces_.size(); functionSpaceNo++)
    {
      Vec3 point = sampledPoints_[samplingPointNo];

      LOG(DEBUG) << "point no " << samplingPointNo << ", point " << point;

      // find the point in the current function space
      element_no_t elementNoLocal = 0;
      int ghostMeshNo = -1;
      Vec2 xi;

      //findPosition(Vec3 point, element_no_t &elementNo, int &ghostMeshNo, std::array<double,D> &xi, bool startSearchInCurrentElement, double xiTolerance)
      bool pointFound = functionSpaces_[functionSpaceNo]->findPosition(point, elementNoLocal, ghostMeshNo, xi, true);


      if (pointFound)
      {
        foundPointNos_.push_back(samplingPointNo);
        elementXis_[samplingPointNo] = std::make_tuple(functionSpaceNo, elementNoLocal, xi);
        LOG(DEBUG) << "point found at el. " << elementNoLocal << ", xi: " << xi;
        break;
      }
      else
      {
        LOG(DEBUG) << "point not found";
      }
    }
  }

  LOG(DEBUG) << "foundPointNos_: " << foundPointNos_;
  LOG(DEBUG) << "elementXis: " << elementXis_;

  // write positions of found sampling points
  std::ofstream fileFoundPoints;
  std::stringstream filenameFoundPoints;

  std::string filenameBase;
  if (filename_.rfind(".") != std::string::npos)
  {
    filenameBase = filename_.substr(0, filename_.rfind("."));
  }
  else
  {
    filenameBase = filename_;
  }

  filenameFoundPoints << filenameBase << "_found." << rankSubset_->ownRankNo() << ".csv";
  Generic::openFile(fileFoundPoints, filenameFoundPoints.str(), true);

  // write found points
  for (int i = 0; i < foundPointNos_.size(); i++)
  {
    int samplingPointNo = foundPointNos_[i];
    Vec3 point = sampledPoints_[samplingPointNo];

    fileFoundPoints << point[0] << ";"
      << point[1] << ";"
      << point[2] << "\n";
  }
  fileFoundPoints.close();

  std::ofstream fileNotFoundPoints;
  std::stringstream filenameNotFoundPoints;


  filenameNotFoundPoints << filenameBase << "_not_found." << rankSubset_->ownRankNo() << ".csv";
  Generic::openFile(fileNotFoundPoints, filenameNotFoundPoints.str(), true);

  // write not found points
  for (int samplingPointNo = 0; samplingPointNo < sampledPoints_.size(); samplingPointNo++)
  {
    if (std::find(foundPointNos_.begin(), foundPointNos_.end(), samplingPointNo) == foundPointNos_.end())
    {
      Vec3 point = sampledPoints_[samplingPointNo];

      fileNotFoundPoints << point[0] << ";"
        << point[1] << ";"
        << point[2] << "\n";
    }
  }

  fileNotFoundPoints.close();
}

template<typename Solver>
void OutputSurface<Solver>::
writeSampledPoints()
{
  if (!ownRankInvolvedInOutput_ || sampledPoints_.empty())
    return;

  LOG(DEBUG) << "writeSampledPoints";

  // get all 2D field variables
  typename DataSurface::FieldVariablesForOutputWriter outputFieldVariables2D
    = this->data_.getFieldVariablesForOutputWriter();

  const int nDofsPerElement = DataSurface::FunctionSpaceFirstFieldVariable::nDofsPerElement();
  const int nComponents = DataSurface::SecondFieldVariable::nComponents();
  if (nComponents != 1)
    LOG(FATAL) << "Output of sampled points is only possible for scalar field variables.";

  std::vector<int> indicesLocal;
  std::vector<double> sampledValuesLocal;
  std::vector<double> sampledGeometryLocal;
  std::vector<int> indicesGlobal;
  std::vector<double> sampledValuesGlobal;
  std::vector<double> sampledGeometryGlobal;

  // read out local point values
  for (int pointNo = 0; pointNo < foundPointNos_.size(); pointNo++)
  {
    int samplingPointNo = foundPointNos_[pointNo];
    int functionSpaceNo = std::get<0>(elementXis_[samplingPointNo]);
    element_no_t elementNoLocal = std::get<1>(elementXis_[samplingPointNo]);
    Vec2 xi = std::get<2>(elementXis_[samplingPointNo]);

    std::shared_ptr<typename DataSurface::SecondFieldVariable> fieldVariable = std::get<2>(outputFieldVariables2D)[functionSpaceNo];

    int nFieldVariables = std::get<2>(outputFieldVariables2D).size();
    if (functionSpaceNo < nFieldVariables)
    {
      // get values in the element
      std::array<double,nDofsPerElement> elementalValues;
      fieldVariable->getElementValues(elementNoLocal, elementalValues);
      LOG(DEBUG) << "get from field variable \"" << fieldVariable->name() << "\", nComponents: " << fieldVariable->nComponents();

      // interpolate value
      double value = fieldVariable->functionSpace()->interpolateValueInElement(elementalValues, xi);
      sampledValuesLocal.push_back(value);

      // get geometry values
      std::array<Vec3,nDofsPerElement> elementalGeometryValues;
      fieldVariable->functionSpace()->geometryField().getElementValues(elementNoLocal, elementalGeometryValues);

      // interpolate value
      Vec3 geometryValue = fieldVariable->functionSpace()->template interpolateValueInElement<3>(elementalGeometryValues, xi);

      for (int componentNo = 0; componentNo < 3; componentNo++)
      {
        sampledGeometryLocal.push_back(geometryValue[componentNo]);
      }

      indicesLocal.push_back(samplingPointNo);
    }
  }

  LOG(DEBUG) << "local values: indicesLocal: " << indicesLocal << ", sampledGeometryLocal: " << sampledGeometryLocal << ", sampledValuesLocal: " << sampledValuesLocal;

  int nRanks = rankSubset_->size();
  int ownRankNo = rankSubset_->ownRankNo();

  // communicate number of found points to rank 0
  std::vector<int> nPointsOnRanks(nRanks);
  int nPointsLocal = foundPointNos_.size();

  MPIUtility::handleReturnValue(MPI_Gather(&nPointsLocal, 1, MPI_INT, nPointsOnRanks.data(), 1, MPI_INT, 0, rankSubset_->mpiCommunicator()), "MPI_Gather");

  // prepare helper fields for gatherv
  std::vector<int> sizesOnRanksIndices(nRanks);
  std::vector<int> sizesOnRanksValues(nRanks);
  std::vector<int> sizesOnRanksGeometry(nRanks);
  std::vector<int> offsetsIndices(nRanks);
  std::vector<int> offsetsValues(nRanks);
  std::vector<int> offsetsGeometry(nRanks);

  int nPointsGlobal = 0;
  if (ownRankNo == 0)
  {
    for (int rankNo = 0; rankNo < nRanks; rankNo++)
    {
      nPointsGlobal += nPointsOnRanks[rankNo];
    }
    // resize global buffers
    indicesGlobal.resize(nPointsGlobal);
    sampledValuesGlobal.resize(nPointsGlobal*nComponents);
    sampledGeometryGlobal.resize(nPointsGlobal*3);

    LOG(DEBUG) << "n global points: " << nPointsGlobal << ", nPointsOnRanks: " << nPointsOnRanks;

    for (int rankNo = 0; rankNo < nRanks; rankNo++)
    {
      sizesOnRanksIndices[rankNo] = nPointsOnRanks[rankNo];
      sizesOnRanksValues[rankNo] = nPointsOnRanks[rankNo] * nComponents;
      sizesOnRanksGeometry[rankNo] = nPointsOnRanks[rankNo] * 3;

      if (rankNo == 0)
      {
        offsetsIndices[rankNo] = 0;
        offsetsValues[rankNo] = 0;
        offsetsGeometry[rankNo] = 0;
      }
      else
      {
        offsetsIndices[rankNo] = offsetsIndices[rankNo-1] + sizesOnRanksIndices[rankNo-1];
        offsetsValues[rankNo] = offsetsValues[rankNo-1] + sizesOnRanksValues[rankNo-1];
        offsetsGeometry[rankNo] = offsetsGeometry[rankNo-1] + sizesOnRanksGeometry[rankNo-1];
      }
    }
  }

  LOG(DEBUG) << "indicesLocal: " << indicesLocal << ", size: " << indicesLocal.size() << "," << indicesGlobal.size() << ", sizesOnRanksValues: "
    << sizesOnRanksValues << ", offsetsValues: " << offsetsValues << ", nPointsLocal*nComponents: " << nPointsLocal*nComponents;

  // communicate actual values
  MPIUtility::handleReturnValue(MPI_Gatherv(indicesLocal.data(), nPointsLocal*nComponents, MPI_DOUBLE,
                                            indicesGlobal.data(), sizesOnRanksValues.data(), offsetsValues.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD), "MPI_Gatherv");
  MPIUtility::handleReturnValue(MPI_Gatherv(sampledValuesLocal.data(), nPointsLocal*nComponents, MPI_DOUBLE,
                                            sampledValuesGlobal.data(), sizesOnRanksValues.data(), offsetsValues.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD), "MPI_Gatherv");
  MPIUtility::handleReturnValue(MPI_Gatherv(sampledGeometryLocal.data(), nPointsLocal*3, MPI_DOUBLE,
                                            sampledGeometryGlobal.data(), sizesOnRanksGeometry.data(), offsetsGeometry.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD), "MPI_Gatherv");

  // on rank 0, write values to file
  if (ownRankNo == 0)
  {
    // sort values
    struct entry_t
    {
      int index;
      Vec3 geometry;
      double value;
    };

    std::vector<entry_t> entries(nPointsGlobal);
    for (int i = 0; i < nPointsGlobal; i++)
    {
      assert(indicesGlobal.size() > i);
      entries[i].index = indicesGlobal[i];

      assert(sampledGeometryGlobal.size() > 3*i+2);
      entries[i].geometry = Vec3{sampledGeometryGlobal[3*i+0], sampledGeometryGlobal[3*i+1], sampledGeometryGlobal[3*i+2]};

      entries[i].value = sampledValuesGlobal[i];
    }

    std::sort(entries.begin(), entries.end(), [](const entry_t &a, const entry_t &b){
      return (a.index == b.index && a.geometry[0] < b.geometry[0]) || (a.index < b.index);
    });

    LOG(DEBUG) << "sorted values: ";
    for (int i = 0; i < nPointsGlobal; i++)
    {
      LOG(DEBUG) << i << ". index=" << entries[i].index << ", geometry: " << entries[i].geometry << ", value: " << entries[i].value;
    }

    //openFile(std::ofstream& file, std::string filename, bool append=false);
    std::ofstream file;
    Generic::openFile(file, filename_, true);

    // write time
    file << currentTime_ << ";";

    // write indices
    for (int i = 0; i < nPointsGlobal; i++)
    {
      file << entries[i].index << ";";
    }

    // write geometry
    for (int i = 0; i < nPointsGlobal; i++)
    {
      file << entries[i].geometry[0] << ";"
        << entries[i].geometry[1] << ";"
        << entries[i].geometry[2] << ";";
    }

    // write values
    for (int i = 0; i < nPointsGlobal; i++)
    {
      file << entries[i].value << ";";
    }
    file.close();

    // at first timestep, write geometry to file
    if (currentTime_ <= 1e-5)
    {
      std::ofstream fileGeometry;
      std::stringstream filenameGeometry;

      std::string filenameBase;
      if (filename_.rfind(".") != std::string::npos)
      {
        filenameBase = filename_.substr(0, filename_.rfind("."));
      }
      else
      {
        filenameBase = filename_;
      }

      filenameGeometry << filenameBase << "_location.csv";
      Generic::openFile(fileGeometry, filenameGeometry.str(), true);

      // write geometry
      for (int i = 0; i < nPointsGlobal; i++)
      {
        fileGeometry << entries[i].geometry[0] << ";"
          << entries[i].geometry[1] << ";"
          << entries[i].geometry[2] << "\n";
      }
      fileGeometry.close();
    }
  }
}

template<typename Solver>
void OutputSurface<Solver>::
run()
{
  initialize();

  solver_.run();

  LOG(DEBUG) << "OutputSurface: writeOutput";
  if (ownRankInvolvedInOutput_)
  {
    outputWriterManager_.writeOutput(data_);
  }
}

template<typename Solver>
void OutputSurface<Solver>::
reset()
{
  solver_.reset();
}

template<typename Solver>
void OutputSurface<Solver>::
setTimeSpan(double startTime, double endTime)
{
  currentTime_ = startTime;
  solver_.setTimeSpan(startTime, endTime);
}

template<typename Solver>
typename OutputSurface<Solver>::Data &OutputSurface<Solver>::
data()
{
  return solver_.data();
}

template<typename Solver>
std::shared_ptr<typename OutputSurface<Solver>::OutputConnectorDataType> OutputSurface<Solver>::
getOutputConnectorData()
{
  return solver_.getOutputConnectorData();
}

}  // namespace OutputWriter
