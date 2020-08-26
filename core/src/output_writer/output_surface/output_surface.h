#pragma once

#include <Python.h>  // has to be the first included header

#include "function_space/function_space.h"
#include "data_management/output_surface/output_surface.h"

namespace OutputWriter
{

/** A class that has own OutputWriters and writes the 2D surface of a 3D mesh.
 */
template<typename Solver>
class OutputSurface :
  public Runnable
{
public:
  typedef typename Solver::FunctionSpace FunctionSpace;
  typedef typename Solver::Data Data;
  typedef typename ::Data::OutputSurface<Data> DataSurface;
  typedef typename Solver::SlotConnectorDataType SlotConnectorDataType;

  //! constructor
  OutputSurface(DihuContext context);

  //! advance simulation by the given time span [startTime_, endTime_]
  void advanceTimeSpan();

  //! initialize time span from specificSettings_
  void initialize();

  //! run solution process
  void run();

  //! reset state
  void reset();

  //! set a new time interval
  void setTimeSpan(double startTime, double endTime);

  //! return the data object of the timestepping scheme
  Data &data();

  //! get the data that will be transferred in the operator splitting to the other term of the splitting
  //! the transfer is done by the slot_connector_data_transfer class
  std::shared_ptr<SlotConnectorDataType> getSlotConnectorData();

protected:

  //! find where the sampling points given in samplingPoints_ are located in the surface mesh
  void initializeSampledPoints();

  //! output the actual sampling points to the file
  void writeSampledPointValues();

  //! write positions of found sampling points
  void writeFoundAndNotFoundPointGeometry();

  DihuContext context_;               //< object that contains the python config for the current context and the global singletons meshManager and solverManager
  Solver solver_;                     //< the contained solver object

  bool initialized_ = false;          //< if this object is initialized
  DataSurface data_;                  //< data object
  bool ownRankInvolvedInOutput_;      //< if the own rank should call the output writer, because surface meshes are output, it can be that the surface is only contained on a subset of ranks
  std::shared_ptr<Partition::RankSubset> rankSubset_;  //< the ranks that are involved in computing the surface
  int timeStepNo_;                    //< time step no for output writer
  double currentTime_;                //< current simulation time for output writer
  double xiTolerance_;                //< tolerance for the xi coordinate used when finding the samplingPoints in the elements of the mesh
  int writeCallCount_ = 0;            //< a counter that increases with the output files

  std::vector<std::shared_ptr<typename ::Data::OutputSurface<Data>::FunctionSpaceFirstFieldVariable>> functionSpaces_;    //< the 2D function spaces for the faces

  struct FoundSampledPoint
  {
    int samplingPointNo;              //< number of the sampling point in the option "sampledPoints"
    Vec3 requestedPosition;           //< the requested position as given in the option "sampledPoints"
    Vec3 position;                    //< the actual position of the point on the mesh
    int functionSpaceNo;              //< no of the function space (which face)
    element_no_t elementNoLocal;      //< element no in the function space that has the point
    Vec2 xi;                          //< xi value within the element where the point is
    double score;                        //< the residual of the point finding algorithm, lower is a better score
  };

  std::map<int,FoundSampledPoint> foundSampledPoints_;    //< all sampled points that have been found locally, key is the samplingPointNo

  std::vector<Vec3> sampledPointsRequestedPositions_;   //< the requested sampling points, as given by the "samplingPoints" option
  std::vector<double> sampledPointsPositionGlobal_;  //< on rank 0, buffer for all geometry of all sampled points
  std::vector<double> valuesGlobal_;  //< on rank 0, buffer for all values
  std::vector<int> partitioningGlobal_; //< on rank 0, buffer for the partitioning values

  std::string filename_;              //< filename of the file to write the sampled points to
  bool updatePointPositions_;         //< if the positions of the sampled points should be updated in every call, this leads to approximately the same global point positions even if the geometry changes
  SeriesWriter seriesWriter_;         //< the series writer object that collects all VTK filenames and creates a collection file that can be loaded by ParaView, for the files that have the EMG values
  SeriesWriter seriesWriterFoundPoints_;      //< the series writer object that collects all VTK filenames and creates a collection file that can be loaded by ParaView, for the files that have the found electrode points
  SeriesWriter seriesWriterNotFoundPoints_;   //< the series writer object that collects all VTK filenames and creates a collection file that can be loaded by ParaView, for the files that have the not found electrode points

  Manager outputWriterManager_;       //< manager object holding all output writers

};

}  // namespace OutputWriter

#include "output_writer/output_surface/output_surface.tpp"
#include "output_writer/output_surface/output_surface_write.tpp"
