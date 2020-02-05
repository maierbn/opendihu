#pragma once

#include <Python.h>  // has to be the first included header

#include <vector>
#include <iostream>
#include <sstream>
#include <list>
#include <functional>
#include <Vc/Allocator>

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
  //! @param approximateExponentialFunction If the exp()-Function should be approximated by the n=1024th series term
  void generateSourceFile(std::string outputFilename, std::string optimizationType, bool approximateExponentialFunction);

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

  //! get additional compile flags that are required to compile the created source file, dependend on the optimizationType, e.g. -fopenmp for "openmp"
  std::string additionalCompileFlags() const;

  //! return the compiler command to use to compile the created source file, e.g. "gcc" or "g++"
  std::string compilerCommand() const;

  //! get the suffix to use for the source file, e.g. ".c" or ".cpp"
  std::string sourceFileSuffix() const;

  struct code_expression_t
  {
    //! type of current syntax node,
    //! type variableName is e.g. "STATES[2]" (split into code="STATES" and arrayIndex=2)
    //! type tree is a list of children, e.g. parantheses: treeChildren = [(type=code, code="("), (type=whatever...), (type=code, code=")"]
    //! type code is just raw code, e.g. operators
    //! type commented_out is not used code that was present in the source file but will not be used for the target file
    enum {variableName, tree, otherCode, commented_out} type;

    std::string code;     //< if type=variableName, code can be STATES, RATES or ALGEBRAIC, if type=code it is normal code
    int arrayIndex;       //< the array index of the state or rate or algebraic

    std::vector<code_expression_t> treeChildren;   //! if type == tree, this contains the children

    //! parse the c++ code line into current struct
    void parse(std::string line);

    //! recursively iterate over every entry with type unequal to tree of the current expression
    void visitLeafs(std::function<void(code_expression_t &expression, bool isFirstVariable)> callback);

    //! implementation of the visitLeafs method
    void visitLeafsCall(std::function<void(CellMLSourceCodeGenerator::code_expression_t &expression, bool isFirstVariable)> callback, bool &isFirstVariable);

    //! visit all nodes, also the tree nodes
    void visitNodes(std::function<void(code_expression_t &expression)> callback);

    //! get a debugging string of the current expression
    std::string getString();

  };

protected:

  //! check if sourceFilename_ is an xml based file and then convert to a c file, updating sourceFilename_
  void convertFromXmlToC();

  //! Scan the given cellml source file and initialize the following:
  //! stateNames_, intermediateNames_, nConstants_ and nIntermediatesInSource_
  //! constantAssignments_, statesInitialValues_,
  //! cellMLCode_
  void parseSourceCodeFile();

  //! write the source file with openmp pragmas in struct-of-array memory ordering
  //! that will be autovectorized by the compiler
  void generateSourceFileSimd(std::string outputFilename);

  //! Write the source file with explicit vectorization using Vc
  //! The file contains the source for only the rhs computation
  void generateSourceFileVc(std::string outputFilename, bool approximateExponentialFunction);

  //! write the source file with explicit vectorization using Vc
  //! The file contains the source for the total solve the rhs computation
  void generateSourceFileSolverExplicitVectorization(std::string outputFilename);

  //! write the source file with openmp support
  void generateSourceFileOpenMP(std::string outputFilename);

  int nInstances_;                             //< number of instances of the CellML problem. Usually it is the number of mesh nodes when a mesh is used. When running in parallel this is the local number of instances without ghosts.
  int nParameters_ = 0;                        //< number of parameters (=CellML name "known") in one instance of the CellML problem

  unsigned int nConstants_ = 0;                //< number of entries in the "CONSTANTS" array
  unsigned int nStates_;                       //< number of states as given in initialize
  unsigned int nIntermediates_;                //< number of intermediates as given in initialize
  unsigned int nIntermediatesInSource_ = 0;    //< number of intermediate values (=CellML name "wanted") in one instance of the CellML problem, as detected from the source file

  std::string compilerCommand_;                //< compiler command that should be used to compile the created source file
  std::string additionalCompileFlags_;         //< additional compile flags that depend on the optimizationType, e.g. -fopenmp for "openmp"
  std::string sourceFileSuffix_;               //< suffix to use for the source file, e.g. ".c" or ".cpp"

  std::vector<int> parametersUsedAsIntermediate_;  //< explicitely defined parameters that will be copied to intermediates, this vector contains the indices of the algebraic array
  std::vector<int> parametersUsedAsConstant_;  //< explicitely defined parameters that will be copied to constants, this vector contains the indices of the constants
  std::vector<double> parameters_;             //< vector of nParameters_*nInstances_ values that will be provided to CellML by the code, given by python config, CellML name: known

  std::vector<std::string> stateNames_;        //< the names for the states as given in the input source file
  std::vector<std::string> intermediateNames_; //< the names for the intermediates as given in the input source file

  std::string sourceFilename_;                 //<file name of provided CellML source file

  std::vector<double> statesInitialValues_;    //< initial values of the states for one instances, as parsed from source file

  std::vector<std::string> constantAssignments_;  //< source code lines where constant variables are assigned

  // contains all the essential parts of the parsed cellml source code
  struct CellMLCode
  {
    std::string header;
    std::vector<code_expression_t> lines;   //< for every line a list of codes
    std::string footer;
  }
  cellMLCode_;
};
