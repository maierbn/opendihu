#pragma once

#include <Python.h>  // has to be the first included header

#include <vector>
#include <iostream>
#include <sstream>
#include <list>

class CellMLSourceCodeGenerator
{
public:
  //! constructor
  CellMLSourceCodeGenerator();

  //! initialize all variables, parses the source code
  void initialize(std::string inputFilename, int nInstances, int nStates, int nIntermediates,
    const std::vector<int> &parametersUsedAsIntermediate, const std::vector<int> &parametersUsedAsConstant,
    const std::vector<double> &parametersInitialValues
  );

  //! generate the source file according to optimizationType
  //! Possible values are: simd vc openmp
  void generateSourceFile(std::string outputFilename, std::string optimizationType) const;

  //! get a reference to the statesInitialValues_ variable
  std::vector<double> &statesInitialValues();

  //! get a reference to the names of the intermediates variables
  const std::vector<std::string> &intermediateNames() const;

  //! get a reference to the names of the state variables
  const std::vector<std::string> &stateNames() const;

  //! get the number of parameters
  const int nParameters() const;

  //! get a reference to the parameters, this allows to change the parameters
  std::vector<double> &parameters();

  //! get the source filename of the initial file (which is inputFilename in initialize)
  const std::string sourceFilename() const;

protected:

  //! Scan the given cellml source file and initialize the following:
  //! stateNames_, intermediateNames_, nConstants_ and nIntermediatesInSource_
  //! constantAssignments_, statesInitialValues_,
  //! cellMLCode_
  void parseSourceCodeFile();

  //! write the source file with openmp pragmas in struct-of-array memory ordering
  //! that will be autovectorized by the compiler
  void generateSourceFileSimdAutovectorization(std::string outputFilename) const;

  //! Write the source file with explicit vectorization using Vc
  //! The file contains the source for only the rhs computation
  void generateSourceFileExplicitVectorization(std::string outputFilename) const;

  //! write the source file with explicit vectorization using Vc
  //! The file contains the source for the total solve the rhs computation
  void generateSourceFileSolverExplicitVectorization(std::string outputFilename) const;

  int nInstances_;                            //< number of instances of the CellML problem. Usually it is the number of mesh nodes when a mesh is used. When running in parallel this is the local number of instances without ghosts.
  int nParameters_ = 0;                       //< number of parameters (=CellML name "known") in one instance of the CellML problem

  int nConstants_ = 0;                        //< number of entries in the "CONSTANTS" array
  int nStates_;                               //< number of states as given in initialize
  int nIntermediates_;                        //< number of intermediates as given in initialize
  int nIntermediatesInSource_ = 0;            //< number of intermediate values (=CellML name "wanted") in one instance of the CellML problem, as detected from the source file

  std::vector<int> parametersUsedAsIntermediate_;  //< explicitely defined parameters that will be copied to intermediates, this vector contains the indices of the algebraic array
  std::vector<int> parametersUsedAsConstant_;  //< explicitely defined parameters that will be copied to constants, this vector contains the indices of the constants
  std::vector<double> parameters_;             //< vector of values that will be provided to CellML by the code, given by python config, CellML name: known

  std::vector<std::string> stateNames_;        //< the names for the states as given in the input source file
  std::vector<std::string> intermediateNames_; //< the names for the intermediates as given in the input source file

  std::string sourceFilename_;                 //<file name of provided CellML source file

  std::vector<double> statesInitialValues_;   //< initial values of the states for one instances, as parsed from source file

  std::vector<std::string> constantAssignments_;  //< source code lines where constant variables are assigned

  struct assignment_t
  {
    std::string code;     //< code can be STATES, RATES or ALGEBRAIC
    int arrayIndex;       //< the array index of the state or rate or algebraic
    enum {variableName, other, commented_out} type;      //< whether it is a variable or other code
  };

  // contains all the essential parts of the parsed cellml source code
  struct CellMLCode
  {
    std::stringstream header;
    std::vector<std::list<assignment_t>> entries;   //< for every line a list of codes
    std::stringstream footer;
  }
  cellMLCode_;
};
