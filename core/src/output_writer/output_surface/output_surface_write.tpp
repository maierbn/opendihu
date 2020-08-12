#include "output_writer/output_surface/output_surface.h"

#include <algorithm>

#include "output_writer/output_surface/output_points.h"

namespace OutputWriter
{

template<typename Solver>
void OutputSurface<Solver>::
writeSampledPointValues()
{
  if (!ownRankInvolvedInOutput_ || sampledPointsRequestedPositions_.empty())
    return;

  LOG(DEBUG) << "writeSampledPointValues";

  // get all 2D field variables
  typename DataSurface::FieldVariablesForOutputWriter outputFieldVariables2D
    = this->data_.getFieldVariablesForOutputWriter();

  const int nDofsPerElement = DataSurface::FunctionSpaceFirstFieldVariable::nDofsPerElement();
  const int nComponents = DataSurface::SecondFieldVariable::nComponents();
  if (nComponents != 1)
    LOG(FATAL) << "Output of sampled points is only possible for scalar field variables.";

  std::vector<int> pointNosLocal;
  std::vector<double> sampledValuesLocal;
  std::vector<double> sampledGeometryLocal;
  std::vector<double> scoresLocal;
  std::vector<int> partitioningLocal;

  std::vector<int> pointNosGlobal;
  std::vector<double> sampledValuesGlobal;
  std::vector<double> sampledGeometryGlobal;
  std::vector<double> scoresGlobal;
  std::vector<int> partitioningGlobal;

  // read out local point values
  // loop over the found points
  for (const std::pair<int,FoundSampledPoint> &pair : foundSampledPoints_)
  {
    int sampledPointNo = pair.first;
    const FoundSampledPoint &foundSampledPoint = pair.second;
    int functionSpaceNo = foundSampledPoint.functionSpaceNo;

    element_no_t elementNoLocal = foundSampledPoint.elementNoLocal;
    Vec2 xi = foundSampledPoint.xi;

    // get the field variable that has the value
    std::shared_ptr<typename DataSurface::SecondFieldVariable> fieldVariable = std::get<2>(outputFieldVariables2D)[functionSpaceNo];

    int nFieldVariables = std::get<2>(outputFieldVariables2D).size();
    if (functionSpaceNo < nFieldVariables)
    {
      // get values in the element
      std::array<double,nDofsPerElement> elementalValues;
      fieldVariable->getElementValues(elementNoLocal, elementalValues);

      LOG(DEBUG) << "get from field variable \"" << fieldVariable->name() << "\", nComponents: " << fieldVariable->nComponents();

      std::shared_ptr<typename ::Data::OutputSurface<Data>::FunctionSpaceFirstFieldVariable> functionSpace = functionSpaces_[functionSpaceNo];

      // interpolate value
      double value = functionSpace->interpolateValueInElement(elementalValues, xi);
      sampledValuesLocal.push_back(value);

      // get geometry values
      std::array<Vec3,nDofsPerElement> elementalGeometryValues;
      functionSpace->geometryField().getElementValues(elementNoLocal, elementalGeometryValues);

      // interpolate value
      Vec3 actualPosition = functionSpace->template interpolateValueInElement<3>(elementalGeometryValues, xi);

      for (int componentNo = 0; componentNo < 3; componentNo++)
      {
        sampledGeometryLocal.push_back(actualPosition[componentNo]);
      }

      pointNosLocal.push_back(sampledPointNo);
      scoresLocal.push_back(foundSampledPoint.score);
    }
  }

  LOG(DEBUG) << "local values: pointNosLocal: " << pointNosLocal << ", sampledGeometryLocal: " << sampledGeometryLocal << ", sampledValuesLocal: " << sampledValuesLocal;

  int nRanks = rankSubset_->size();
  int ownRankNo = rankSubset_->ownRankNo();

  // communicate number of found points to rank 0
  std::vector<int> nPointsOnRanks(nRanks);
  int nPointsLocal = pointNosLocal.size();
  partitioningLocal.resize(nPointsLocal, ownRankNo);

  MPIUtility::handleReturnValue(MPI_Gather(&nPointsLocal, 1, MPI_INT, nPointsOnRanks.data(), 1, MPI_INT, 0, rankSubset_->mpiCommunicator()), "MPI_Gather");

  // prepare helper fields for gatherv
  static std::vector<int> sizesOnRanksIndices;
  static std::vector<int> sizesOnRanksValues;
  static std::vector<int> sizesOnRanksGeometry;
  static std::vector<int> offsetsIndices;
  static std::vector<int> offsetsValues;
  static std::vector<int> offsetsGeometry;

  sizesOnRanksIndices.resize(nRanks);
  sizesOnRanksValues.resize(nRanks);
  sizesOnRanksGeometry.resize(nRanks);
  offsetsIndices.resize(nRanks);
  offsetsValues.resize(nRanks);
  offsetsGeometry.resize(nRanks);

  int nPointsGlobal = 0;
  if (ownRankNo == 0)
  {
    for (int rankNo = 0; rankNo < nRanks; rankNo++)
    {
      nPointsGlobal += nPointsOnRanks[rankNo];
    }
    // resize global buffers
    pointNosGlobal.resize(nPointsGlobal);
    sampledValuesGlobal.resize(nPointsGlobal*nComponents, 0);
    sampledGeometryGlobal.resize(nPointsGlobal*3, 0);
    scoresGlobal.resize(nPointsGlobal, 0);
    partitioningGlobal.resize(nPointsGlobal, 0);

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

  LOG(DEBUG) << "pointNosLocal: " << pointNosLocal << ", size: " << pointNosLocal.size() << "," << pointNosGlobal.size() << ", sizesOnRanksValues: "
    << sizesOnRanksValues << ", offsetsValues: " << offsetsValues << ", nPointsLocal*nComponents: " << nPointsLocal*nComponents;

  // communicate actual values
  MPIUtility::handleReturnValue(MPI_Gatherv(pointNosLocal.data(), nPointsLocal*nComponents, MPI_INT,
                                            pointNosGlobal.data(), sizesOnRanksValues.data(), offsetsValues.data(), MPI_INT, 0, rankSubset_->mpiCommunicator()), "MPI_Gatherv");
  MPIUtility::handleReturnValue(MPI_Gatherv(sampledValuesLocal.data(), nPointsLocal*nComponents, MPI_DOUBLE,
                                            sampledValuesGlobal.data(), sizesOnRanksValues.data(), offsetsValues.data(), MPI_DOUBLE, 0, rankSubset_->mpiCommunicator()), "MPI_Gatherv");
  MPIUtility::handleReturnValue(MPI_Gatherv(sampledGeometryLocal.data(), nPointsLocal*3, MPI_DOUBLE,
                                            sampledGeometryGlobal.data(), sizesOnRanksGeometry.data(), offsetsGeometry.data(), MPI_DOUBLE, 0, rankSubset_->mpiCommunicator()), "MPI_Gatherv");
  MPIUtility::handleReturnValue(MPI_Gatherv(scoresLocal.data(), nPointsLocal, MPI_DOUBLE,
                                            scoresGlobal.data(), sizesOnRanksValues.data(), offsetsValues.data(), MPI_DOUBLE, 0, rankSubset_->mpiCommunicator()), "MPI_Gatherv");
  MPIUtility::handleReturnValue(MPI_Gatherv(partitioningLocal.data(), nPointsLocal, MPI_INT,
                                            partitioningGlobal.data(), sizesOnRanksValues.data(), offsetsValues.data(), MPI_INT, 0, rankSubset_->mpiCommunicator()), "MPI_Gatherv");

  // on rank 0, write values to file
  if (ownRankNo == 0)
  {
    // sort values
    struct entry_t
    {
      int pointNo;
      Vec3 geometry;
      double value;
      double score;
      int partition;
    };

    std::vector<entry_t> entries(nPointsGlobal);
    for (int i = 0; i < nPointsGlobal; i++)
    {
      assert(pointNosGlobal.size() > i);
      entries[i].pointNo = pointNosGlobal[i];

      assert(sampledGeometryGlobal.size() > 3*i+2);
      entries[i].geometry = Vec3{sampledGeometryGlobal[3*i+0], sampledGeometryGlobal[3*i+1], sampledGeometryGlobal[3*i+2]};

      entries[i].value = sampledValuesGlobal[i];
      entries[i].score = scoresGlobal[i];
      entries[i].partition = partitioningGlobal[i];
    }

    // sort collected points from all ranks by their sampledPointNo, then by score
    std::sort(entries.begin(), entries.end(), [](const entry_t &a, const entry_t &b)
    {
      return (a.pointNo < b.pointNo)
        || (a.pointNo == b.pointNo && a.score < b.score)
        || (a.pointNo == b.pointNo && a.score == b.score && a.geometry[0] < b.geometry[0]);
    });

#ifndef NDEBUG
    LOG(DEBUG) << "sorted values: ";
    for (int i = 0; i < nPointsGlobal; i++)
    {
      LOG(DEBUG) << i << ". pointNo: " << entries[i].pointNo << ", geometry: " << entries[i].geometry
        << ", value: " << entries[i].value << ", score: " << entries[i].score;
    }
#endif

    // remove duplicates (same pointNo) and store in fields
    sampledPointsPositionGlobal_.clear();
    valuesGlobal_.clear();
    partitioningGlobal_.clear();
    sampledPointsPositionGlobal_.reserve(nPointsGlobal*3);
    valuesGlobal_.reserve(nPointsGlobal);
    partitioningGlobal_.reserve(nPointsGlobal);

    int lastI = -1;
    for (int i = 0; i < nPointsGlobal; i++)
    {
      if (lastI >= 0)
        if (entries[i].pointNo == entries[lastI].pointNo)
        {
          LOG(DEBUG) << "skip point " << i << ", " << entries[i].pointNo << " at " << entries[i].geometry << ", score " << entries[i].score
            << ", last was " << lastI << ", " << entries[lastI].pointNo << " at " << entries[lastI].geometry << ", score " << entries[lastI].score;
          continue;
        }

      sampledPointsPositionGlobal_.push_back(entries[i].geometry[0]);
      sampledPointsPositionGlobal_.push_back(entries[i].geometry[1]);
      sampledPointsPositionGlobal_.push_back(entries[i].geometry[2]);
      valuesGlobal_.push_back(entries[i].value);
      partitioningGlobal_.push_back(entries[i].partition);
      lastI = i;
    }

    // create a filename for the vtp file
    std::stringstream vtpFile;
    if (filename_.rfind(".") != std::string::npos)
    {
      vtpFile << filename_.substr(0, filename_.rfind("."));
    }
    else
    {
      vtpFile << filename_;
    }
    vtpFile << "_" << std::setw(7) << std::setfill('0') << writeCallCount_ << ".vtp";
    writeCallCount_++;

    // write result to files
    OutputPoints::writeCsvFile(filename_, currentTime_, sampledPointsPositionGlobal_, valuesGlobal_);
    OutputPoints::writeVtpFile(vtpFile.str(), currentTime_, sampledPointsPositionGlobal_, valuesGlobal_, 1, partitioningGlobal_, "EMG");

    seriesWriter_.registerNewFile(vtpFile.str(), currentTime_);
  }
}

template<typename Solver>
void OutputSurface<Solver>::
writeFoundAndNotFoundPointGeometry()
{
  // write positions of found sampling points
  std::string filenameBase;
  if (filename_.rfind(".") != std::string::npos)
  {
    filenameBase = filename_.substr(0, filename_.rfind("."));
  }
  else
  {
    filenameBase = filename_;
  }

  // write found points of current rank
  std::stringstream filenameFoundPointsCsv;
  std::stringstream filenameFoundPointsVtp;
  filenameFoundPointsCsv << filenameBase << "_found." << rankSubset_->ownRankNo() << ".csv";
  filenameFoundPointsVtp << filenameBase << "_found." << rankSubset_->ownRankNo() << ".vtp";

  std::vector<double> foundPointsGeometry;
  std::vector<double> empty;

  // loop over the found points
  for (const std::pair<int,FoundSampledPoint> &pair : foundSampledPoints_)
  {
    Vec3 point = pair.second.requestedPosition;
    foundPointsGeometry.push_back(point[0]);
    foundPointsGeometry.push_back(point[1]);
    foundPointsGeometry.push_back(point[2]);
  }

  int ownRankNo = DihuContext::ownRankNoCommWorld();
  std::vector<int> partitioning(foundPointsGeometry.size(), ownRankNo);

  OutputPoints::writeCsvFile(filenameFoundPointsCsv.str(), currentTime_, foundPointsGeometry, empty);
  OutputPoints::writeVtpFile(filenameFoundPointsVtp.str(), currentTime_, foundPointsGeometry, empty, 1, partitioning, "EMG");

  // write not found points of current rank
  std::stringstream filenameNotFoundPointsCsv;
  std::stringstream filenameNotFoundPointsVtp;
  filenameNotFoundPointsCsv << filenameBase << "_not_found." << rankSubset_->ownRankNo() << ".csv";
  filenameNotFoundPointsVtp << filenameBase << "_not_found." << rankSubset_->ownRankNo() << ".vtp";

  std::vector<double> notFoundPointsGeometry;

  for (int samplingPointNo = 0; samplingPointNo < sampledPointsRequestedPositions_.size(); samplingPointNo++)
  {
    if (foundSampledPoints_.find(samplingPointNo) == foundSampledPoints_.end())
    {
      Vec3 point = sampledPointsRequestedPositions_[samplingPointNo];
      notFoundPointsGeometry.push_back(point[0]);
      notFoundPointsGeometry.push_back(point[1]);
      notFoundPointsGeometry.push_back(point[2]);
    }
  }

  partitioning.resize(notFoundPointsGeometry.size(), ownRankNo);
  OutputPoints::writeCsvFile(filenameNotFoundPointsCsv.str(), currentTime_, notFoundPointsGeometry, empty);
  OutputPoints::writeVtpFile(filenameNotFoundPointsVtp.str(), currentTime_, notFoundPointsGeometry, empty, 1, partitioning, "EMG");

  seriesWriterFoundPoints_.registerNewFile(filenameFoundPointsVtp.str(), currentTime_);
  seriesWriterNotFoundPoints_.registerNewFile(filenameNotFoundPointsVtp.str(), currentTime_);
}

}  // namespace OutputWriter
