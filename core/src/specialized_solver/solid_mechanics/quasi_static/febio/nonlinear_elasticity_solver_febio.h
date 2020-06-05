#pragma once

#include <Python.h>  // has to be the first included header

#include "data_management/specialized_solver/quasi_static_nonlinear_elasticity_febio.h"
#include "data_management/output_connector_data.h"
#include "output_writer/manager.h"

namespace TimeSteppingScheme
{

/** A specialized solver for 3D linear elasticity, incompressible Mooney-Rivlin.
 *  This class simply computes the static problem. The activation for a muscle material is ignored in this class.
 * 
 *  This is also the base class for QuasiStaticNonlinearElasticitySolverFebio. The 
 *  QuasiStaticNonlinearElasticitySolverFebio class computes the quasi static problem using the febio muscle material in a timestepping scheme.
 *
 *  This class produces febio files that can only be computed by febio3
 */
class NonlinearElasticitySolverFebio :
  public Runnable
{
public:
  typedef FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>, BasisFunction::LagrangeOfOrder<1>> FunctionSpace;
  typedef ::Data::QuasiStaticNonlinearElasticityFebio Data;
  typedef FieldVariable::FieldVariable<FunctionSpace,1> FieldVariableType;
  typedef Data::OutputConnectorDataType OutputConnectorDataType;

  //! constructor
  NonlinearElasticitySolverFebio(DihuContext context, std::string solverName="NonlinearElasticitySolverFebio");

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
  //! the transfer is done by the output_connector_data_transfer class
  std::shared_ptr<OutputConnectorDataType> getOutputConnectorData();

  //! output the given data for debugging
  std::string getString(std::shared_ptr<OutputConnectorDataType> data);

  //! determine if febio2 was installed and is available on the command line
  bool isFebioAvailable();

protected:

  //! create the febio_input.feb file which contains the problem for febio to solve
  virtual void createFebioInputFile();

  //! load the file that was created by the febio simulation
  void loadFebioOutputFile();

  //! run the febio program on the generated input file
  void runFebio();

  //! communicate elemtal values such as global node nos and activation values to rank 0
  void communicateElementValues(std::vector<double> &activationValuesGlobal, std::vector<int> &nodeNosGlobal);

  //! communicate all nodal values to rank 0 to be written to the febio input file
  void communicateNodeValues(std::vector<double> &nodePositionValuesGlobal);

  DihuContext context_;    //< object that contains the python config for the current context and the global singletons meshManager and solverManager

  OutputWriter::Manager outputWriterManager_; //< manager object holding all output writer
  Data data_;                               //< data object

  std::string durationLogKey_;              //< key with with the duration of the computation is written to the performance measurement log
  std::string solverName_;                  //< the name of the config, i.e. "NonlinearElasticitySolverFebio" or "QuasiStaticNonlinearElasticitySolverFebio"
  std::string problemDescription_;          //< text to appear in the febio input file
  std::string problemTitle_;                //< title to appear in the febio input file

  PythonConfig specificSettings_;           //< python object containing the value of the python config dict with corresponding key
  double activationFactor_;                 //< factor with which to multiply activation
  std::vector<double> materialParameters_;  //< the material parameters c0, c1 and k for Mooney-Rivlin material
  std::vector<element_no_t> tractionElementNos_;     //< element nos of elements that will have traction
  Vec3 tractionVector_;                     //< the traction vector to apply on the traction surfaces
  
  double endTime_;     //< end time of current time step
  bool initialized_;   //< if initialize() was already called
};

}  // namespace

