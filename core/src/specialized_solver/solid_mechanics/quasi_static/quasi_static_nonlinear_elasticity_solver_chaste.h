#pragma once

#include <Python.h>  // has to be the first included header

#include "data_management/specialized_solver/quasi_static_nonlinear_elasticity_chaste.h"

namespace TimeSteppingScheme
{

/** This solver is an adapter to the chaste software framework and uses the hyperelasticity implementation of chaste.
 *  It solves the nonlinear finite elasticity problem with Mooney-Rivlin material, for either 2D or 3D.
 *  Because Chaste is not able to solve nonlinear elasticity in parallel, nor solve anything else than the quasi-static case,
 *  integration in opendihu is not complete. This class should be deleted now.
 *
 */
template<int D>
class QuasiStaticNonlinearElasticitySolverChaste :
  public Runnable
{
public:
  typedef FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<D>, BasisFunction::LagrangeOfOrder<2>> FunctionSpaceType;
  typedef Data::QuasiStaticNonlinearElasticityChaste<FunctionSpaceType> Data;
  typedef FieldVariable::FieldVariable<FunctionSpaceType,1> FieldVariableType;
  typedef ::Data::SlotConnectorData<FunctionSpaceType,1> SlotConnectorDataType;

  //! constructor
  QuasiStaticNonlinearElasticitySolverChaste(DihuContext context);

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
  //! the transfer is done by the slot_connector_data_transfer class
  std::shared_ptr<SlotConnectorDataType> getSlotConnectorData();

  //! output the given data for debugging
  std::string getString(std::shared_ptr<SlotConnectorDataType> data);

protected:

  DihuContext context_;    //< object that contains the python config for the current context and the global singletons meshManager and solverManager

  OutputWriter::Manager outputWriterManager_; //< manager object holding all output writer
  Data data_;                 //< data object

  std::string durationLogKey_;   //< key with with the duration of the computation is written to the performance measurement log
  std::shared_ptr<FunctionSpaceType> functionSpace_;  //< the function space with quadratic Lagrange basis functions, as created in the opendihu code. This will be passed on to chaste.

  bool initialized_;   //< if this object was already initialized
  PythonConfig specificSettings_;    //< python object containing the value of the python config dict with corresponding key
  double endTime_;     //< end time of current time step
  double maximumActiveStress_;    //< parameter value of the maximum active stress, this is the scaling factor of the activation value to get the active stress tensor
  double strainScalingCurveWidth_;   //< width of a parabola that scales the stress dependend on the relative sarcomere length
};

}  // namespace

#include "specialized_solver/solid_mechanics/quasi_static/quasi_static_nonlinear_elasticity_solver_chaste.tpp"
