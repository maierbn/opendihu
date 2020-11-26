#pragma once

#include <Python.h>  // has to be the first included header
#include <iostream>
#include <vector>
#include <atomic>

#include "control/types.h"
#include "output_writer/generic.h"

namespace OutputWriter
{

struct BoundingBox
{
  Vec3 min, max;   // minimum,maximum value of bounding box
  BoundingBox();
};

struct MegaMolWriterContext
{
  std::vector<Vec3> geometryFieldVectorValues;    //< value buffer for geometry values, stored as Vec3
  std::vector<double> geometryFieldValues;        //< value buffer for geometry values, stored continuously

  std::array<double,6> globalBoundingBoxValues;   //< value buffer for global bounding box values
  std::array<double,6> localBoundingBoxValues;    //< value buffer for local bounding box values

  std::vector<double> vmValues;                   //< value buffer of Vm values, for both EMG and fibers
  std::vector<double> emgValues;                  //< value buffer for emg values, only for 3D data set
  std::vector<double> transmembraneFlowValues;    //< value buffer for transmembraneFlowValues, only for 3D data set
  std::vector<int> partitioning;                  //< value buffer for
  std::array<int,3> nPointsPerCoordinateDirection;
  std::vector<int> rankNo;                        //< value buffer for a field that contains the own global rank no for all ranks

  double approximateDistanceBetweenFibers;        //< estimated distance between fiber
  int nodeOffsetOnOwnComputeNode;                 //< buffer for the offset value

  //! clear all buffers and reset values
  void clear();
};

class MegaMol : public Generic
{
public:

  //! constructor
  MegaMol(DihuContext context, PythonConfig specificSettings, std::shared_ptr<Partition::RankSubset> rankSubset = nullptr);

  //! destructor
  //~MegaMol();

  //! write out solution to given filename, if timeStepNo is not -1, this value will be part of the filename
  template<typename DataType>
  void write(DataType &data, int timeStepNo = -1, double currentTime = -1, int callCountIncrement = 1);

#ifdef HAVE_ADIOS
  struct adios_writer_t
  {
    MegaMolWriterContext megaMolWriterContext;  //< all variables that are combined and collected during the calls until combineNInstances_ calls have been made

    std::shared_ptr<adios2::IO> adiosIo;        //< io object that is associated with defined variables, creates the engine
    std::shared_ptr<adios2::Engine> engine;     //< adios writer
    std::atomic_int nOpenWriters;               //< this writer could be in multiple instances, called at the same time. Then just the first open should really open the file and the last close should close it. This variables counts the number of opens.
  };

private:

#if defined(HAVE_MEGAMOL)
  void notifyMegaMol();
#endif

  //! write all variables with ADIOS
  void writeAdiosVariables();

  //! write the geometry information
  void writeGeometryField();

  //! write the local and global bounding boxes
  void writeBoundaryBoxes();

  //! write the payload data consisting of scalar fields
  void writeScalarFields();

  static std::map<std::string,std::array<std::shared_ptr<adios_writer_t>, 2>> adiosWriters_;   //< the writers for front and back buffer, key in the map is the filenameBase_
  static std::map<std::string,std::shared_ptr<adios_writer_t>> adiosWriter_;   //< the currently used writer, key in the map is the filenameBase_
  int currentOpenWriterIndex_ = 0;   //< which writer the last opened is, 0 or 1

  std::shared_ptr<adios2::Variable<double>> adiosFieldVariableGeometry_;            //< the adios field variable for the geometry, with name "xyz"
  std::shared_ptr<adios2::Variable<double>> adiosFieldVariableVm_;                  //< the adios scalar field variable with name "vm"
  std::shared_ptr<adios2::Variable<double>> adiosFieldVariableEmg_;                 //< the adios scalar field variable with name "emg"
  std::shared_ptr<adios2::Variable<double>> adiosFieldVariableTransmembraneFlow_;   //< the adios scalar field variable with name "transmembraneFlow"
  std::shared_ptr<adios2::Variable<int>> adiosNPointsPerCoordinateDirection_;       //< the adios scalar field variable, with name "nPointsPerCoordinateDirection"
  std::shared_ptr<adios2::Variable<int>> offsetVariable_;                           //< the adios scalar field variable, with name "offset"
  std::shared_ptr<adios2::Variable<double>> localBoundingBoxVariable_;              //< the adios variable containing the local bounding box information, with name "localBoundingBox"
  std::shared_ptr<adios2::Variable<double>> globalBoundingBoxVariable_;             //< the adios variable containing the global bounding box information, with name "globalBoundingBox"
  std::shared_ptr<adios2::Variable<double>> globalRadiusVariable_;                  //< the adios variable for the radius for visualization of spheres, with name "global_radius"
  std::shared_ptr<adios2::Variable<int>> globalNumberOfNodesVariable_;              //< the adios variable for the global number of nodes
  std::shared_ptr<adios2::Variable<int>> rankNoVariable_;                           //< the adios variable for the rank no, for every node in the field

  double globalRadius_;         //< the global radius to be used for the visualization
  int nNodesGlobal_ = 0;        //< the global number of nodes or the particle count

  std::string currentFilename_; //< the file to which is currently being written, the file is not yet ready
  std::string lastFilename_;    //< the last used filename for output

  int combineNInstances_;       //< number of calls to MegaMol::write that should be combined into one timestep, this is needed when the output writer is given in a MultipleInstances environment
  bool useFrontBackBuffer_;     //< if the two buffers, *_0 and *_1 should be used, otherwise is uses a standard file naming scheme with increasing counter number
#endif

};

} // namespace

#include "output_writer/megamol/megamol.tpp"
