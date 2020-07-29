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
  void writeSampledPoints();

  DihuContext context_;               //< object that contains the python config for the current context and the global singletons meshManager and solverManager
  Solver solver_;                     //< the contained solver object

  bool initialized_ = false;          //< if this object is initialized
  DataSurface data_;                  //< data object
  bool ownRankInvolvedInOutput_;      //< if the own rank should call the output writer, because surface meshes are output, it can be that the surface is only contained on a subset of ranks
  std::shared_ptr<Partition::RankSubset> rankSubset_;  //< the ranks that are involved in computing the surface
  int timeStepNo_;                    //< time step no for output writer
  double currentTime_;                //< current simulation time for output writer

  std::vector<std::shared_ptr<typename ::Data::OutputSurface<Data>::FunctionSpaceFirstFieldVariable>> functionSpaces_;    //< the 2D function spaces for the faces
  std::vector<Vec3> sampledPoints_;   //< if set, the point coordinates where to sample the 2D surface mesh
  std::vector<int> foundPointNos_;    //< the point nos of the found points
  std::map<int,std::tuple<int,element_no_t,Vec2>> elementXis_;   //< the function space no (which face), element no and xi value for each sampling point on the local domain
  std::string filename_;              //< filename of the file to write the sampled points to

  Manager outputWriterManager_;       //< manager object holding all output writers

};

}  // namespace OutputWriter

#include "output_writer/output_surface/output_surface.tpp"
