#include "output_writer/output_surface/output_surface.h"

#include <algorithm>

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
                                            pointNosGlobal.data(), sizesOnRanksValues.data(), offsetsValues.data(), MPI_INT, 0, MPI_COMM_WORLD), "MPI_Gatherv");
  MPIUtility::handleReturnValue(MPI_Gatherv(sampledValuesLocal.data(), nPointsLocal*nComponents, MPI_DOUBLE,
                                            sampledValuesGlobal.data(), sizesOnRanksValues.data(), offsetsValues.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD), "MPI_Gatherv");
  MPIUtility::handleReturnValue(MPI_Gatherv(sampledGeometryLocal.data(), nPointsLocal*3, MPI_DOUBLE,
                                            sampledGeometryGlobal.data(), sizesOnRanksGeometry.data(), offsetsGeometry.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD), "MPI_Gatherv");
  MPIUtility::handleReturnValue(MPI_Gatherv(scoresLocal.data(), nPointsLocal, MPI_DOUBLE,
                                            scoresGlobal.data(), sizesOnRanksValues.data(), offsetsValues.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD), "MPI_Gatherv");
  MPIUtility::handleReturnValue(MPI_Gatherv(partitioningLocal.data(), nPointsLocal, MPI_INT,
                                            partitioningGlobal.data(), sizesOnRanksValues.data(), offsetsValues.data(), MPI_INT, 0, MPI_COMM_WORLD), "MPI_Gatherv");

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
    sampledPointsPositionGlobal_.reserve(nPointsGlobal);
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

      sampledPointsPositionGlobal_.push_back(entries[i].geometry);
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
    writeCsvFile(filename_, currentTime_, sampledPointsPositionGlobal_, valuesGlobal_);
    writeVtpFile(vtpFile.str(), currentTime_, sampledPointsPositionGlobal_, valuesGlobal_, partitioningGlobal_);

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

  std::vector<Vec3> foundPointsGeometry;
  std::vector<double> empty;

  // loop over the found points
  for (const std::pair<int,FoundSampledPoint> &pair : foundSampledPoints_)
  {
    Vec3 point = pair.second.requestedPosition;
    foundPointsGeometry.push_back(point);
  }

  int ownRankNo = DihuContext::ownRankNoCommWorld();
  std::vector<int> partitioning(foundPointsGeometry.size(), ownRankNo);

  writeCsvFile(filenameFoundPointsCsv.str(), currentTime_, foundPointsGeometry, empty);
  writeVtpFile(filenameFoundPointsVtp.str(), currentTime_, foundPointsGeometry, empty, partitioning);

  // write not found points of current rank
  std::stringstream filenameNotFoundPointsCsv;
  std::stringstream filenameNotFoundPointsVtp;
  filenameNotFoundPointsCsv << filenameBase << "_not_found." << rankSubset_->ownRankNo() << ".csv";
  filenameNotFoundPointsVtp << filenameBase << "_not_found." << rankSubset_->ownRankNo() << ".vtp";

  std::vector<Vec3> notFoundPointsGeometry;

  for (int samplingPointNo = 0; samplingPointNo < sampledPointsRequestedPositions_.size(); samplingPointNo++)
  {
    if (foundSampledPoints_.find(samplingPointNo) == foundSampledPoints_.end())
    {
      Vec3 point = sampledPointsRequestedPositions_[samplingPointNo];
      notFoundPointsGeometry.push_back(point);
    }
  }

  partitioning.resize(notFoundPointsGeometry.size(), ownRankNo);
  writeCsvFile(filenameNotFoundPointsCsv.str(), currentTime_, notFoundPointsGeometry, empty);
  writeVtpFile(filenameNotFoundPointsVtp.str(), currentTime_, notFoundPointsGeometry, empty, partitioning);

  seriesWriterFoundPoints_.registerNewFile(filenameFoundPointsVtp.str(), currentTime_);
  seriesWriterNotFoundPoints_.registerNewFile(filenameNotFoundPointsVtp.str(), currentTime_);
}

template<typename Solver>
void OutputSurface<Solver>::
writeCsvFile(std::string filename, double currentTime, const std::vector<Vec3> &geometry, const std::vector<double> &values)
{
  std::ofstream file;
  Generic::openFile(file, filename, true);  // append to file

  int nPointsGlobal = geometry.size();

  // in first timestep, write header
  if (currentTime <= 1e-5)
  {
    file << "#timestamp;t;n_points";
    for (int pointNo = 0; pointNo < nPointsGlobal; pointNo++)
    {
      file << ";p" << pointNo << "_x;p" << pointNo << "_y;p" << pointNo << "_z";
    }

    if (!values.empty())
    {
      for (int pointNo = 0; pointNo < nPointsGlobal; pointNo++)
      {
        file << ";p" << pointNo << "_value";
      }
    }
    file << std::endl;
  }

  // write timestamp and time
  // time stamp
  auto t = std::time(nullptr);
  auto tm = *std::localtime(&t);
  file << StringUtility::timeToString(&tm) << ";"
    << currentTime_ << ";" << nPointsGlobal;

  // write geometry
  for (int i = 0; i < nPointsGlobal; i++)
  {
    file << ";" << geometry[i][0]
      << ";" << geometry[i][1]
      << ";" << geometry[i][2];
  }

  // write values
  if (!values.empty())
  {
    for (int i = 0; i < nPointsGlobal; i++)
    {
      file << ";" << values[i];
    }
  }
  file << std::endl;
  file.close();
}

template<typename Solver>
void OutputSurface<Solver>::
writeVtpFile(std::string filename, double currentTime, const std::vector<Vec3> &geometry, const std::vector<double> &values,
             const std::vector<int> &partitioning)
{
  int nPointsGlobal = geometry.size();

  // transform values vector to vector with contiguous entries
  static std::vector<double> geometryValues;
  geometryValues.resize(nPointsGlobal*3);
  for (int pointNo = 0; pointNo < nPointsGlobal; pointNo++)
  {
    geometryValues[3*pointNo + 0] = geometry[pointNo][0];
    geometryValues[3*pointNo + 1] = geometry[pointNo][1];
    geometryValues[3*pointNo + 2] = geometry[pointNo][2];
  }

  std::ofstream file;
  Generic::openFile(file, filename, false);  // do not append to existing file

  // transform current time to string
  std::vector<double> time(1, this->currentTime_);
  std::string stringTime = Paraview::encodeBase64Float(time.begin(), time.end());

  file << "<?xml version=\"1.0\"?>" << std::endl
    << "<!-- " << DihuContext::versionText() << " " << DihuContext::metaText()
    << ", currentTime: " << this->currentTime_ << " -->" << std::endl
    << "<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\">" << std::endl    // intel cpus are LittleEndian
    << std::string(1, '\t') << "<PolyData>" << std::endl
    << std::string(2, '\t') << "<FieldData>" << std::endl
    << std::string(3, '\t') << "<DataArray type=\"Float32\" Name=\"Time\" NumberOfTuples=\"1\" format=\"binary\" >" << std::endl
    << std::string(4, '\t') << stringTime << std::endl
    << std::string(3, '\t') << "</DataArray>" << std::endl
    << std::string(2, '\t') << "</FieldData>" << std::endl;

  file << std::string(2, '\t') << "<Piece NumberOfPoints=\"" << nPointsGlobal << "\" NumberOfVerts=\"0\" "
    << "NumberOfLines=\"0\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">" << std::endl
    << std::string(3, '\t') << "<PointData>" << std::endl;


  // write emg values
  if (!values.empty())
  {
    file << std::string(4, '\t') << "<DataArray "
        << "Name=\"EMG\" "
        << "type=\"Float32\" "
        << "NumberOfComponents=\"1\" "
        << "format=\"binary\" >" << std::endl << std::string(5, '\t')
        << Paraview::encodeBase64Float(values.begin(), values.end())
        << std::endl << std::string(4, '\t') << "</DataArray>" << std::endl;
  }

  // write partitioning values
  file << std::string(4, '\t') << "<DataArray "
      << "Name=\"partitioning\" "
      << "type=\"Int32\" "
      << "NumberOfComponents=\"1\" "
      << "format=\"binary\" >" << std::endl << std::string(5, '\t')
      << Paraview::encodeBase64Int32(partitioning.begin(), partitioning.end())
      << std::endl << std::string(4, '\t') << "</DataArray>" << std::endl;

  file << std::string(3, '\t') << "</PointData>" << std::endl
    << std::string(3, '\t') << "<CellData>" << std::endl
    << std::string(3, '\t') << "</CellData>" << std::endl
    << std::string(3, '\t') << "<Points>" << std::endl
    << std::string(4, '\t') << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"binary"
    << "\" >" << std::endl << std::string(5, '\t');

  // write geometry values
  file << Paraview::encodeBase64Float(geometryValues.begin(), geometryValues.end());

  file << std::endl << std::string(4, '\t') << "</DataArray>" << std::endl
    << std::string(3, '\t') << "</Points>" << std::endl
    << std::string(3, '\t') << "<Verts></Verts>" << std::endl
    << std::string(3, '\t') << "<Lines></Lines>" << std::endl
    << std::string(3, '\t') << "<Strips></Strips>" << std::endl
    << std::string(3, '\t') << "<Polys></Polys>" << std::endl
    << std::string(2, '\t') << "</Piece>" << std::endl
    << std::string(1, '\t') << "</PolyData>" << std::endl
    << "</VTKFile>" << std::endl;

  file.close();
}

}  // namespace OutputWriter
