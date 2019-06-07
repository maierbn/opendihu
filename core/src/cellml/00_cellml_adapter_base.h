#pragma once

#include <Python.h>  // has to be the first included header

#include <vector>

#include "control/dihu_context.h"
#include "output_writer/manager.h"
#include "function_space/function_space.h"

/** This is the base class of the CellmlAdapter, that handles common functionality.
 * nStates: number of states in one instance of the CellML problem
 * 
 *  Naming:
 *   Intermediate (opendihu) = KNOWN (OpenCMISS) = Algebraic (OpenCOR)
 *   Parameter (opendihu, OpenCMISS) = KNOWN (OpenCMISS), in OpenCOR also algebraic
 *   Constant - these are constants that are only present in the source files
 *   State: state variable
 *   Rate: the time derivative of the state variable, i.e. the increment value in an explicit Euler stepping
 */
template <int nStates, int nIntermediates_, typename FunctionSpaceType>
class CellmlAdapterBase
{
public:

  typedef FieldVariable::FieldVariable<FunctionSpaceType,nIntermediates_> FieldVariableTypeIntermediates;

  //! constructor from context
  CellmlAdapterBase(DihuContext context, bool noNewOutputWriter);

  //! constructor from context
  CellmlAdapterBase(DihuContext context);

  //! destructor
  ~CellmlAdapterBase();

  //! return the compile-time constant number of state variables of one instance that will be integrated
  static constexpr int nComponents();

  //! load model, use settings given in context
  void initialize();

  //! create generic function space with given number of instances
  void initializeFromNInstances(int nInstances);
  
  //! set initial values as given in python config
  template<typename FunctionSpaceType2>
  bool setInitialValues(std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType2,nStates>> initialValues);

  //! return false because the object is independent of mesh type
  bool knowsMeshType();

  //! return the mesh
  std::shared_ptr<FunctionSpaceType> functionSpace();

  //! get number of instances, number of intermediates and number of parameters
  void getNumbers(int &nInstances, int &nIntermediates, int &nParameters);

  //! get a vector with the names of the states
  void getStateNames(std::vector<std::string> &stateNames);

  //! get the outputStateIndex value, which is the index of the state that should be used further in an operator splitting scheme, for electrophysiology application this is the states of Vm
  int outputStateIndex();

  //! get the outputIntermediateIndex, which is the index of the intermediate that should be used for further computation in the operator splitting
  int outputIntermediateIndex();

  //! get the prefactor value, i.e. the factor with which the solution will be scaled before the transfer in an operator splitting scheme
  double prefactor();

  //! get the const number of intermediates
  constexpr int nIntermediates() const;

  //! return vector of the intermediate values
  std::shared_ptr<FieldVariableTypeIntermediates> intermediates();

protected:

  //! scan the given cellml source file for initial values that are given by dummy assignments (OpenCMISS) or directly (OpenCOR). This also sets nParameters_, nConstants_ and nIntermediatesFromSource_
  virtual bool scanSourceFile(std::string sourceFilename, std::array<double,nStates> &statesInitialValues) = 0;

  DihuContext context_;    ///< object that contains the python config for the current context and the global singletons meshManager and solverManager
  PythonConfig specificSettings_;    ///< python object containing the value of the python config dict with corresponding key
  OutputWriter::Manager outputWriterManager_; ///< manager object holding all output writer

  std::shared_ptr<FunctionSpaceType> functionSpace_;    ///< a mesh, there are as many instances of the same CellML problem as there are nodes in the mesh

  int nInstances_;         ///< number of instances of the CellML problem. Usually it is the number of mesh nodes when a mesh is used. When running in parallel this is the local number of instances without ghosts.
  int nParameters_ = 0;    ///< number of parameters (=CellML name "known") in one instance of the CellML problem
  int nIntermediatesInSource_ = 0; ///< number of intermediate values (=CellML name "wanted") in one instance of the CellML problem, as detected from the source file
  int nConstants_ = 0;     ///< number of entries in the "CONSTANTS" array
   
  int outputStateIndex_ = 0;   ///< the index of the state that should be used further in an operator splitting scheme, for electrophysiology application this is the states of Vm
  int outputIntermediateIndex_ = 0;  ///< the index of the intermediates that should be transferred to further solvers, analogous to outputStateIndex
  double prefactor_ = 0;       ///< the factor with which the solution will be scaled before the transfer in an operator splitting scheme
  int internalTimeStepNo_ = 0; ///< the counter how often the right hand side was called

  //std::vector<double> states_;    ///< vector of states, that are computed by rhsRoutine, this is not needed as member variable, because the states are directly stored in the Petsc Vecs of the solving time stepping scheme
  //std::vector<double> rates_;     ///< vector of rates, that are computed by rhsRoutine, this is not needed as member variable, because the states are directly stored in the Petsc Vecs of the solving time stepping scheme
  std::vector<double> parameters_; ///< vector of values that will be provided to CellML by the code, given by python config, CellML name: known
  //std::vector<double> intermediates_;    ///< vector of intermediate values in DAE system. These can be computed directly from the actual states at any time. Gets computed by rhsRoutine from states, together with rates. OpenCMISS name is intermediate, CellML name: wanted
  std::array<double,nStates> statesInitialValues_;  ///< initial values of the states for one instances, as parsed from source file
  
  std::shared_ptr<FieldVariableTypeIntermediates> intermediates_;   ///< field variable that hold the intermediate values

  std::vector<int> parametersUsedAsIntermediate_;  ///< explicitely defined parameters that will be copied to intermediates, this vector contains the indices of the algebraic array
  std::vector<int> parametersUsedAsConstant_;  ///< explicitely defined parameters that will be copied to constants, this vector contains the indices of the constants 
  
  std::vector<std::string> stateNames_;    ///< the specifier for the states as given in the input source file
  
  std::string sourceFilename_; ///<file name of provided CellML right hand side routine
  bool inputFileTypeOpenCMISS_;   ///< if the input file that is being parsed is from OpenCMISS and not from OpenCOR
};

#include "cellml/00_cellml_adapter_base.tpp"
