#pragma once

#include <Python.h>  // has to be the first included header

#include <vector>
#include <iostream>
#include <sstream>
#include <list>
#include <functional>
#include <vc_or_std_simd.h>

#ifndef HAVE_STDSIMD      // only if we are using Vc, it is not necessary for std::simd
#include <Vc/Allocator>
#endif

class CellmlSourceCodeGeneratorBase
{
public:
  //! constructor
  CellmlSourceCodeGeneratorBase();

  //! initialize the algebraicNames, stateNames and constantNames by parsing the source code. There will be an error if the number of states/algebraics does not match
  void initializeNames(std::string inputFilename, int nInstances, int nStates, int nAlgebraics);

  //! initialize all variables, parses the source code
  void initializeSourceCode(
    const std::vector<int> &parametersUsedAsAlgebraic, const std::vector<int> &parametersUsedAsConstant,
    std::vector<double> &parametersInitialValues, int maximumNumberOfParameters, double *parameterValues
  );

  //! generate the source file according to optimizationType
  //! Possible values are: simd vc openmp
  //! @param approximateExponentialFunction If the exp()-Function should be approximated by the n=1024th series term
  void generateSourceFile(std::string outputFilename, std::string optimizationType, bool approximateExponentialFunction);

  //! get a reference to the statesInitialValues_ variable
  std::vector<double> &statesInitialValues();

  //! get a reference to the names of the algebraics variables
  const std::vector<std::string> &algebraicNames() const;

  //! get a reference to the names of the state variables
  const std::vector<std::string> &stateNames() const;

  //! get a reference to the names of the constants
  const std::vector<std::string> &constantNames() const;

  //! get the number of parameters
  const int nParameters() const;

  //! get the source filename of the initial file (which is inputFilename in initialize)
  const std::string sourceFilename() const;

  //! get additional compile flags that are required to compile the created source file, dependend on the optimizationType, e.g. -fopenmp for "openmp"
  std::string additionalCompileFlags() const;

  //! return the compiler command to use to compile the created source file, e.g. "gcc" or "g++"
  std::string compilerCommand() const;

  //! get the suffix to use for the source file, e.g. ".c" or ".cpp"
  std::string sourceFileSuffix() const;

protected:

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
    void visitLeafsCall(std::function<void(CellmlSourceCodeGeneratorBase::code_expression_t &expression, bool isFirstVariable)> callback, bool &isFirstVariable);

    //! visit all nodes, also the tree nodes
    void visitNodes(std::function<void(code_expression_t &expression)> callback);

    //! get a debugging string of the current expression
    std::string getString();

  };

  //! check if sourceFilename_ is an xml based file and then convert to a c file, updating sourceFilename_
  void convertFromXmlToC();

  //! Scan the cellml source file and initialize the following:
  //! stateNames_, algebraicNames_, nConstants_ and nAlgebraicsInSource_
  void parseNamesInSourceCodeFile();

  //! Scan the given cellml source file and initialize the following:
  //! constantAssignments_, statesInitialValues_,
  //! cellMLCode_
  void parseSourceCodeFile();

  //! Generate the rhs code for a single instance. This is needed for computing the equilibrium of the states.
  void generateSingleInstanceCode();



  int nInstances_;                             //< number of instances of the CellML problem. Usually it is the number of mesh nodes when a mesh is used. When running in parallel this is the local number of instances without ghosts.
  int nParameters_ = 0;                        //< number of parameters (=CellML name "known") in one instance of the CellML problem

  unsigned int nConstants_ = 0;                //< number of entries in the "CONSTANTS" array
  unsigned int nStates_;                       //< number of states as given in initialize
  unsigned int nAlgebraics_;                //< number of algebraics as given in initialize
  unsigned int nAlgebraicsInSource_ = 0;    //< number of algebraic values (=CellML name "wanted") in one instance of the CellML problem, as detected from the source file

  std::string compilerCommand_;                //< compiler command that should be used to compile the created source file
  std::string additionalCompileFlags_;         //< additional compile flags that depend on the optimizationType, e.g. -fopenmp for "openmp"
  std::string sourceFileSuffix_;               //< suffix to use for the source file, e.g. ".c" or ".cpp"

  std::vector<int> parametersUsedAsAlgebraic_;  //< explicitely defined parameters that will be copied to algebraics, this vector contains the indices of the algebraic array
  std::vector<int> parametersUsedAsConstant_;  //< explicitely defined parameters that will be copied to constants, this vector contains the indices of the constants

  std::vector<std::string> stateNames_;        //< the names for the states as given in the input source file
  std::vector<std::string> algebraicNames_; //< the names for the algebraics as given in the input source file
  std::vector<std::string> constantNames_;     //< the names of the constants

  std::string sourceFilename_;                 //< file name of provided CellML source file
  std::string singleInstanceCode_;             //< c++ code that computes the rhs for a single instance, as in the original source file. This is needed to compute the equilibrium of states.

  std::vector<double> statesInitialValues_;    //< initial values of the states for one instance, as parsed from source file

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
