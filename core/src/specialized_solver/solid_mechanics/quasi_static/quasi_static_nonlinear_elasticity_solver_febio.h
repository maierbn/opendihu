#pragma once

#include <Python.h>  // has to be the first included header

#include "data_management/specialized_solver/quasi_static_nonlinear_elasticity_febio.h"
#include "specialized_solver/solid_mechanics/quasi_static/output_connector_data_type.h"
#include "output_writer/manager.h"

namespace TimeSteppingScheme
{

/** A specialized solver for 3D linear elasticity, as quasi-static timestepping scheme (a new static solution every timestep)
  */
class QuasiStaticNonlinearElasticitySolverFebio :
  public Runnable
{
public:
  typedef FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>, BasisFunction::LagrangeOfOrder<2>> FunctionSpace;
  typedef ::Data::QuasiStaticNonlinearElasticityFebio Data;
  typedef FieldVariable::FieldVariable<FunctionSpace,1> FieldVariableType;
  typedef ElasticitySolverOutputConnectorDataType<FieldVariableType> OutputConnectorDataType;

  //! constructor
  QuasiStaticNonlinearElasticitySolverFebio(DihuContext context);

  //! advance simulation by the given time span, data in solution is used, afterwards new data is in solution
  void advanceTimeSpan();

  //! initialize components of the simulation
  void initialize();

  //! dummy method, set endTime as current output time
  void setTimeSpan(double startTime, double endTime);

  //! run the simulation
  void run();

  //! reset state
  void reset();

  //! return the data object
  Data &data();

  //! get the data that will be transferred in the operator splitting to the other term of the splitting
  //! the transfer is done by the solution_vector_mapping class
  OutputConnectorDataType getOutputConnectorData();

  //! output the given data for debugging
  std::string getString(OutputConnectorDataType &data);

protected:

  //! create the febio_input.feb file which contains the problem for febio to solve
  void createFebioInputFile();

  //! load the file that was created by the febio simulation
  void loadFebioOutputFile();

  //! run the febio program on the generated input file
  void runFebio();

  DihuContext context_;    ///< object that contains the python config for the current context and the global singletons meshManager and solverManager

  OutputWriter::Manager outputWriterManager_; ///< manager object holding all output writer
  Data data_;                 ///< data object

  std::string durationLogKey_;   ///< key with with the duration of the computation is written to the performance measurement log

  PythonConfig specificSettings_;    ///< python object containing the value of the python config dict with corresponding key
  double endTime_;     ///< end time of current time step
  bool initialized_;   ///< if initialize() was already called
};

}  // namespace

