#pragma once

#include <Python.h>  // has to be the first included header

#include "interfaces/runnable.h"
#include "function_space/function_space.h"
#include "basis_function/basis_function.h"
#include "data_management/specialized_solver/quasi_static_hyperelasticity.h"
#include "output_writer/manager.h"

namespace TimeSteppingScheme
{

/** This solver is for the nonlinear finite elasticity problem with Mooney-Rivlin material in 3D.
  */
class QuasiStaticHyperelasticitySolver :
  public Runnable
{
public:
  typedef FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>, BasisFunction::LagrangeOfOrder<2>> DisplacementsFunctionSpace;
  typedef FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>, BasisFunction::LagrangeOfOrder<1>> PressureFunctionSpace;
  typedef Data::QuasiStaticHyperelasticity<PressureFunctionSpace,DisplacementsFunctionSpace> Data;

  typedef FieldVariable::FieldVariable<DisplacementsFunctionSpace,3> DisplacementsFieldVariableType;
  typedef FieldVariable::FieldVariable<PressureFunctionSpace,1> PressureFieldVariableType;
  typedef FieldVariable::FieldVariable<DisplacementsFunctionSpace,6> StressFieldVariableType;

  //! constructor
  QuasiStaticHyperelasticitySolver(DihuContext context);

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

  //! return whether the underlying discretizableInTime object has a specified mesh type and is not independent of the mesh type
  bool knowsMeshType();

  //! return the data object
  Data &data();

protected:

  DihuContext context_;    ///< object that contains the python config for the current context and the global singletons meshManager and solverManager

  OutputWriter::Manager outputWriterManager_; ///< manager object holding all output writer
  Data data_;                 ///< data object

  std::string durationLogKey_;   ///< key with with the duration of the computation is written to the performance measurement log
  std::shared_ptr<DisplacementsFunctionSpace> displacementsFunctionSpace_;  ///< the function space with quadratic Lagrange basis functions, used for discretization of displacements
  std::shared_ptr<PressureFunctionSpace> pressureFunctionSpace_;  ///< the function space with linear Lagrange basis functions, used for discretization of pressure

  bool initialized_;   ///< if this object was already initialized
  PythonConfig specificSettings_;    ///< python object containing the value of the python config dict with corresponding key
  double endTime_;     ///< end time of current time step
  double c0_;    ///< first Mooney-Rivlin parameter
  double c1_;   ///< second Mooney-Rivlin parameter
};

}  // namespace
