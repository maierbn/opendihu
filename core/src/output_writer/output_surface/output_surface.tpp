#include "output_writer/output_surface/output_surface.h"

#include <algorithm>

namespace OutputWriter
{

template<typename Solver>
OutputSurface<Solver>::
OutputSurface(DihuContext context) :
  context_(context["OutputSurface"]), solver_(context_),
  data_(context_), ownRankInvolvedInOutput_(true), timeStepNo_(0), currentTime_(0.0), updatePointPositions_(false)
{

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
  sampledPointsRequestedPositions_ = PythonUtility::convertFromPython<std::vector<Vec3>>::get(samplingPointsPy);

  if (!sampledPointsRequestedPositions_.empty())
  {
    filename_ = specificSettings.getOptionString("filename", "out/sampledPoints.csv");
    updatePointPositions_ = specificSettings.getOptionBool("updatePointPositions", false);
    xiTolerance_ = specificSettings.getOptionDouble("xiTolerance", 0.3, PythonUtility::NonNegative);
    enableCsvFile_ = specificSettings.getOptionBool("enableCsvFile", true);
    enableVtpFile_ = specificSettings.getOptionBool("enableVtpFile", true);
    enableGeometryInCsvFile_ = specificSettings.getOptionBool("enableGeometryInCsvFile", true);
  }

  LOG(DEBUG) << "OutputSurface: initialize output writers";

  // initialize output writer to use smaller rank subset that only contains the ranks that have parts of the surface
  // if the last argument is not given, by default the common rank subset would be used
  if (ownRankInvolvedInOutput_)
  {
    rankSubset_ = data_.functionSpace()->meshPartition()->rankSubset();
    outputWriterManager_.initialize(context_, specificSettings, rankSubset_);
    std::ofstream file;
    Generic::openFile(file, filename_);  // recreate and truncate file
    file.close();

    initializeSampledPoints();

    // write positions of found sampling points
    if (!sampledPointsRequestedPositions_.empty())
      writeFoundAndNotFoundPointGeometry();
  }

  initialized_ = true;
}

template<typename Solver>
void OutputSurface<Solver>::
initializeSampledPoints()
{
  if (sampledPointsRequestedPositions_.empty())
    return;

  LOG(DEBUG) << "initialize sampled points";
  foundSampledPoints_.clear();

  // get the 2D function spaces
  this->data_.getFunctionSpaces(functionSpaces_);

  for (int functionSpaceNo = 0; functionSpaceNo < functionSpaces_.size(); functionSpaceNo++)
  {
    LOG(DEBUG) << "functionSpace " << functionSpaceNo << "/" << functionSpaces_.size() << ": " << functionSpaces_[functionSpaceNo]->meshName();
  }


  const int nDofsPerElement = DataSurface::FunctionSpaceFirstFieldVariable::nDofsPerElement();
  // now we have a 2D function space for each face in functionSpaces

  // loop over sampling points and find them in the function spaces
  for (int samplingPointNo = 0; samplingPointNo < sampledPointsRequestedPositions_.size(); samplingPointNo++)
  {
    Vec3 point = sampledPointsRequestedPositions_[samplingPointNo];

    // loop over function spaces for the faces
    for (int functionSpaceNo = 0; functionSpaceNo < functionSpaces_.size(); functionSpaceNo++)
    {
      std::shared_ptr<typename ::Data::OutputSurface<Data>::FunctionSpaceFirstFieldVariable> functionSpace = functionSpaces_[functionSpaceNo];

      LOG(DEBUG) << "point no " << samplingPointNo << ", point " << point;

      // find the point in the current function space
      element_no_t elementNoLocal = 0;
      int ghostMeshNo = -1;
      Vec2 xi;

      double residual;
      bool searchedAllElements;
      //findPosition(Vec3 point, element_no_t &elementNo, int &ghostMeshNo, std::array<double,D> &xi, bool startSearchInCurrentElement, double &residual, bool &searchedAllElements, double xiTolerance)
      bool pointFound = functionSpace->findPosition(point, elementNoLocal, ghostMeshNo, xi, true, residual, searchedAllElements, xiTolerance_);

      if (pointFound)
      {
        // fill data structure with information about found point
        FoundSampledPoint foundSampledPoint;
        foundSampledPoint.samplingPointNo = samplingPointNo;
        foundSampledPoint.requestedPosition = sampledPointsRequestedPositions_[samplingPointNo];
        foundSampledPoint.functionSpaceNo = functionSpaceNo;
        foundSampledPoint.elementNoLocal = elementNoLocal;
        foundSampledPoint.xi = xi;

        // determine actual position on the mesh
        std::array<Vec3,nDofsPerElement> elementalGeometryValues;
        functionSpace->geometryField().getElementValues(elementNoLocal, elementalGeometryValues);
        foundSampledPoint.position = functionSpace->template interpolateValueInElement<3>(elementalGeometryValues, xi);

        //foundSampledPoint.requestedPosition = foundSampledPoint.position ;

        // compute score, smaller is better
        foundSampledPoint.score = MathUtility::distance<3>(foundSampledPoint.position, foundSampledPoint.requestedPosition);
        //foundSampledPoint.score = residual;

        // store it
        foundSampledPoints_[samplingPointNo] = foundSampledPoint;
        LOG(DEBUG) << "point found at el. " << elementNoLocal << ", xi: " << xi << ", score: " << foundSampledPoint.score;
        break;
      }
      else
      {
        LOG(DEBUG) << "point not found";
      }
    }
  }
}

template<typename Solver>
void OutputSurface<Solver>::
advanceTimeSpan()
{
  solver_.advanceTimeSpan();

  LOG(DEBUG) << "OutputSurface: writeOutput, ownRankInvolvedInOutput_: " << ownRankInvolvedInOutput_;

  // if the own rank has a part of the surface that will be written
  if (ownRankInvolvedInOutput_)
  {
    // write 2D surface files
    outputWriterManager_.writeOutput(data_, timeStepNo_++);

    // write out values at points
    writeSampledPointValues();

    // write positions of found sampling points
    if (updatePointPositions_)
    {
      initializeSampledPoints();
      writeFoundAndNotFoundPointGeometry();
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
std::shared_ptr<typename OutputSurface<Solver>::SlotConnectorDataType> OutputSurface<Solver>::
getSlotConnectorData()
{
  return solver_.getSlotConnectorData();
}

}  // namespace OutputWriter
