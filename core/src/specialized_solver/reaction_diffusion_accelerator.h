#pragma once

#include <Python.h>  // has to be the first included header
#include "time_stepping_scheme/time_stepping_scheme.h"
#include "interfaces/runnable.h"
#include "control/dihu_context.h"

namespace TimeSteppingScheme
{

/** 
  */
template <int nStates, int nIntermediates=9>
class ReactionDiffusionAccelerator :
  public CellMLAdapter<nStates,nIntermediates,FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>,BasisFunction::LagrangeOfOrder<1>>>, 
  public TimeSteppingScheme,
  public Runnable
{
public:
  typedef FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>,BasisFunction::LagrangeOfOrder<1>> FunctionSpaceType;
  typedef std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>> TransferableSolutionDataType;
  typedef typename Data::ReactionDiffusionAccelerator<FunctionSpaceType,nStates> Data;

  //! constructor
  ReactionDiffusionAccelerator(DihuContext context);

  //! advance simulation by the given time span, data in solution is used, afterwards new data is in solution
  void advanceTimeSpan();

  //! initialize components of the simulation
  void initialize();
  
  //! setting initial values from CellML code
  bool setInitialValues(std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>> initialValues);

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
  Data data_;  ///< the data object of the multidomain solver which stores all field variables and matrices

  OutputWriter::Manager outputWriterManager_; ///< manager object holding all output writer

  std::shared_ptr<Partition::RankSubset> rankSubset_;  ///< the rankSubset for all involved ranks

  std::string durationLogKey_;   ///< key with with the duration of the computation is written to the performance measurement log

  bool initialized_;   ///< if this object was already initialized
  PythonConfig specificSettings_;    ///< python object containing the value of the python config dict with corresponding key
};

}  // namespace

#include "specialized_solver/reaction_diffusion_accelerator.tpp"
