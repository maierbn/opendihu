#pragma once

#include <Python.h>  // has to be the first included header

#include "interfaces/runnable.h"
#include "time_stepping_scheme/00_time_stepping_scheme.h"
#include "data_management/specialized_solver/muscle_contraction_solver.h"
#include "specialized_solver/solid_mechanics/dynamic_hyperelasticity/dynamic_hyperelasticity_solver.h"

/** Solve the incompressible, transversely isotropic Mooney-Rivlin material with active stress contribution.
 */
class MuscleContractionSolver :
  public Runnable,
  public TimeSteppingScheme::TimeSteppingScheme
{
public:
  typedef DynamicHyperelasticitySolver<Equation::SolidMechanics::MooneyRivlinIncompressible3D> DynamicHyperelasticitySolverType;

  //! make the DisplacementsFunctionSpace of the DynamicHyperelasticitySolver class available
  typedef typename DynamicHyperelasticitySolverType::DisplacementsFunctionSpace FunctionSpace;

  //! define the type of the data object,
  typedef typename DynamicHyperelasticitySolverType::Data Data;

  //! Define the type of data that will be transferred between solvers when there is a coupling scheme.
  //! Usually you define this type in the "Data" class and reuse it here.
  typedef typename Data::OutputConnectorDataType OutputConnectorDataType;

  //! constructor, gets the DihuContext object which contains all python settings
  MuscleContractionSolver(DihuContext context);

  //! advance simulation by the given time span [startTime_, endTime_] (set by setTimeSpan(), take a look at time_stepping_scheme/00_time_stepping_scheme.h)
  void advanceTimeSpan();

  //! initialize time span from specificSettings_
  void initialize();

  //! run solution process, this first calls initialize() and then run()
  void run();

  //! reset state of this object, such that a new initialize() is necessary ("uninitialize")
  void reset();

  //! return the data object of the timestepping scheme, with the call to this method the output writers get the data to create their output files
  Data &data();

  //! Get the data that will be transferred in the operator splitting or coupling to the other term of the splitting/coupling.
  //! the transfer is done by the output_connector_data_transfer class
  std::shared_ptr<OutputConnectorDataType> getOutputConnectorData();

protected:

  DynamicHyperelasticitySolverType dynamicHyperelasticitySolver_;   //< the dynamic hyperelasticity solver that solves for the contraction
};
