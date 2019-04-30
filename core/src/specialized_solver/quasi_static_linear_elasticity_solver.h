#pragma once

#include <Python.h>  // has to be the first included header
#include "time_stepping_scheme/time_stepping_scheme_ode.h"
#include "interfaces/runnable.h"
#include "data_management/finite_element_method/finite_elements.h"
#include "control/dihu_context.h"
#include "partition/rank_subset.h"

namespace TimeSteppingScheme
{

/** A specialized solver for linear elasticity, as quasi-static timestepping scheme
  */
template<typename FiniteElementMethod>
class QuasiStaticLinearElasticitySolver :
  public Runnable
{
public:
  typedef typename FiniteElementMethod::FunctionSpace FunctionSpace;
  typedef typename Data::FiniteElements<typename FiniteElementMethodDiffusion::FunctionSpace,3,Equation::Static::LinearElasticity> DataLinearElasticityType;
  typedef Data::QuasiStaticLinearElasticitySolver<DataLinearElasticityType> Data;
  typedef typename FiniteElementsData::FieldVariableType FieldVariableType;
  typedef std::shared_ptr<FieldVariableType> TransferableSolutionDataType;

  //! constructor
  QuasiStaticLinearElasticitySolver(DihuContext context);

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

  //! get the data that will be transferred in the operator splitting to the other term of the splitting
  //! the transfer is done by the solution_vector_mapping class
  TransferableSolutionDataType getSolutionForTransfer();

  //! output the given data for debugging
  std::string getString(TransferableSolutionDataType &data);

protected:

  DihuContext context_;    ///< object that contains the python config for the current context and the global singletons meshManager and solverManager

  OutputWriter::Manager outputWriterManager_; ///< manager object holding all output writer
  Data data_;                 ///< data object

  FiniteElementMethod finiteElementMethodLinearElasticity_;   ///< the finite element object that solves the linear elasticity equation

  std::string durationLogKey_;   ///< key with with the duration of the computation is written to the performance measurement log

  bool initialized_;   ///< if this object was already initialized
  PythonConfig specificSettings_;    ///< python object containing the value of the python config dict with corresponding key
  double endTime_;     ///< end time of current time step
};

}  // namespace

#include "specialized_solver/quasi_static_linear_elasticity.tpp"
