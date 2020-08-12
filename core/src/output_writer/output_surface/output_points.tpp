#include "output_writer/output_surface/output_points.h"

#include <algorithm>

namespace OutputWriter
{

template<int nComponents>
void OutputPoints::
outputPointsVtp(std::string filename, double currentTime, const std::vector<Vec3> &geometry,
                const std::vector<VecD<nComponents>> &vectorValues, std::string fieldVariableName,
                std::shared_ptr<Partition::RankSubset> rankSubset)
{
  // communicate number of points to rank 0
  int nPointsLocal = geometry.size();

  // transform values vector to vector with contiguous entries, for geometry
  geometryValuesLocal_.resize(geometry.size() * 3);
  for (int pointNo = 0; pointNo < geometry.size(); pointNo++)
  {
    for (int componentNo = 0; componentNo < 3; componentNo++)
    {
      geometryValuesLocal_[3*pointNo + componentNo] = geometry[pointNo][componentNo];
    }
  }

  // transform values vector to vector with contiguous entries
  valuesLocal_.resize(vectorValues.size() * nComponents);
  for (int pointNo = 0; pointNo < vectorValues.size(); pointNo++)
  {
    for (int componentNo = 0; componentNo < nComponents; componentNo++)
    {
      valuesLocal_[nComponents*pointNo + componentNo] = vectorValues[pointNo][componentNo];
    }
  }


  int nRanks = rankSubset->size();
  int ownRankNo = rankSubset->ownRankNo();

  // fill local partitioning vector with ownRankNo entries
  partitioningLocal_.resize(nPointsLocal, ownRankNo);

  nPointsOnRanks_.resize(nRanks);
  MPIUtility::handleReturnValue(MPI_Gather(&nPointsLocal, 1, MPI_INT, nPointsOnRanks_.data(), 1, MPI_INT, 0, rankSubset->mpiCommunicator()), "MPI_Gather");

  // prepare helper fields for gatherv
  sizesOnRanksValues_.resize(nRanks);
  sizesOnRanksGeometry_.resize(nRanks);
  sizesOnRanksPartitioning_.resize(nRanks);

  offsetsValues_.resize(nRanks);
  offsetsGeometry_.resize(nRanks);
  offsetsPartitioning_.resize(nRanks);

  int nPointsGlobal = 0;
  if (ownRankNo == 0)
  {
    for (int rankNo = 0; rankNo < nRanks; rankNo++)
    {
      nPointsGlobal += nPointsOnRanks_[rankNo];
    }
    // resize global buffers
    valuesGlobal_.resize(nPointsGlobal*nComponents, 0);
    geometryValuesGlobal_.resize(nPointsGlobal*3, 0);
    partitioningGlobal_.resize(nPointsGlobal, 0);

    LOG(DEBUG) << "n global points: " << nPointsGlobal << ", nPointsOnRanks_: " << nPointsOnRanks_;

    for (int rankNo = 0; rankNo < nRanks; rankNo++)
    {
      sizesOnRanksValues_[rankNo] = nPointsOnRanks_[rankNo] * nComponents;
      sizesOnRanksGeometry_[rankNo] = nPointsOnRanks_[rankNo] * 3;
      sizesOnRanksPartitioning_[rankNo] = nPointsOnRanks_[rankNo];

      if (rankNo == 0)
      {
        offsetsValues_[rankNo] = 0;
        offsetsGeometry_[rankNo] = 0;
        offsetsPartitioning_[rankNo] = 0;
      }
      else
      {
        offsetsValues_[rankNo] = offsetsValues_[rankNo-1] + sizesOnRanksValues_[rankNo-1];
        offsetsGeometry_[rankNo] = offsetsGeometry_[rankNo-1] + sizesOnRanksGeometry_[rankNo-1];
        offsetsPartitioning_[rankNo] = offsetsPartitioning_[rankNo-1] + sizesOnRanksPartitioning_[rankNo-1];
      }
    }
  }

  // communicate actual values
  MPIUtility::handleReturnValue(MPI_Gatherv(valuesLocal_.data(), nPointsLocal*nComponents, MPI_DOUBLE,
                                            valuesGlobal_.data(), sizesOnRanksValues_.data(), offsetsValues_.data(), MPI_DOUBLE, 0, rankSubset->mpiCommunicator()), "MPI_Gatherv");
  MPIUtility::handleReturnValue(MPI_Gatherv(geometryValuesLocal_.data(), nPointsLocal*3, MPI_DOUBLE,
                                            geometryValuesGlobal_.data(), sizesOnRanksGeometry_.data(), offsetsGeometry_.data(), MPI_DOUBLE, 0, rankSubset->mpiCommunicator()), "MPI_Gatherv");
  MPIUtility::handleReturnValue(MPI_Gatherv(partitioningLocal_.data(), nPointsLocal, MPI_INT,
                                            partitioningGlobal_.data(), sizesOnRanksPartitioning_.data(), offsetsPartitioning_.data(), MPI_INT, 0, rankSubset->mpiCommunicator()), "MPI_Gatherv");

  // on rank 0, write values to file
  if (ownRankNo == 0)
  {
    if (filename.find(".vtp") == std::string::npos)
    {
      filename = filename + std::string(".vtp");
    }

    // write result to files
    writeVtpFile(filename, currentTime, geometryValuesGlobal_, valuesGlobal_, nComponents, partitioningGlobal_, fieldVariableName);
  }
}

}  // namespace OutputWriter
