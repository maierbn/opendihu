#pragma once

#include <Python.h>  // has to be the first included header

#include "control/types.h"
#include "partition/rank_subset.h"

namespace OutputWriter
{

/** This class can write *.vtp or *.csv files with points, a single field variable and partitioning information.
 *  There are static functions that write the files serially, or the parallel function that first gathers the values from all ranks
 *  and then writes the file on rank 0.
 *
 *  This class is used by the OutputSurface writer that outputs values for EMG electrodes and by the DirichletBoundaryConditions that outputs Dirichlet BC values.
 */
class OutputPoints
{
public:

  //! write a VTK file with geometry and optionlly values, this has to be called in parallel, it communicates all given values to rank 0 which writes the file
  template<int nComponents>
  void outputPointsVtp(std::string filename, double currentTime, const std::vector<Vec3> &geometry,
                       const std::vector<VecD<nComponents>> &values, std::string fieldVariableName,
                       std::shared_ptr<Partition::RankSubset> rankSubset);

  //! write a VTK file with geometry and optionally values, this is completely serial
  //! @param values This contains the contiguous array of struct representation of the values with nComponents components. If empty, no values are written, only geometry
  //! @param geometry geometry values in array of struct representation with 3 components
  //! @param nComponents the number of components of the values in values
  static void writeVtpFile(std::string filename, double currentTime, const std::vector<double> &geometry,
                           const std::vector<double> &values, int nComponents, const std::vector<int> &partitioning, std::string fieldVariableName);

  //! write a csv file with geometry and optionally values
  static void writeCsvFile(std::string filename, double currentTime, const std::vector<double> &geometry, const std::vector<double> &values, bool writeGeometry=true);

protected:

  std::vector<int> nPointsOnRanks_;           //< how many points there are on every rank

  std::vector<double> valuesLocal_;           //< buffer for local value
  std::vector<double> valuesGlobal_;          //< buffer for global values
  std::vector<int> sizesOnRanksValues_;       //< field needed for MPI_Gatherv, the number of values on every rank
  std::vector<int> offsetsValues_;            //< field needed for MPI_Gatherv, offsets of the arrays on every rank

  std::vector<double> geometryValuesLocal_;   //< buffer for local values of geometry
  std::vector<double> geometryValuesGlobal_;  //< buffer for global values of geometry
  std::vector<int> sizesOnRanksGeometry_;     //< field needed for MPI_Gatherv, the number of values on every rank
  std::vector<int> offsetsGeometry_;          //< field needed for MPI_Gatherv, offsets of the arrays on every rank

  std::vector<int> partitioningLocal_;        //<  buffer for local values of partition no
  std::vector<int> partitioningGlobal_;    //<  buffer for global values of partition no
  std::vector<int> sizesOnRanksPartitioning_; //< field needed for MPI_Gatherv, the number of values on every rank
  std::vector<int> offsetsPartitioning_;      //< field needed for MPI_Gatherv, offsets of the arrays on every rank
};

}  // namespace OutputWriter

#include "output_writer/output_surface/output_points.tpp"
