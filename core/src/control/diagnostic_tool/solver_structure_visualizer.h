#pragma once

#include <Python.h>  // has to be the first included header

#include <iostream>
#include <vector>
#include <memory>
#include <functional>
#include <tuple>

//! forward declarations
namespace Data{
template<typename FunctionSpaceType, int nComponents1, int nComponents2>
class OutputConnectorData;
}
class OutputConnection;

/** Class that collects information about all nested solvers and produces a diagram as txt file.
 *
 */
class SolverStructureVisualizer
{
public:

  //! constructor
  SolverStructureVisualizer();

  //! add a solver to the diagram
  void addSolver(std::string name, bool hasInternalConnectionToFirstNestedSolver=false);

  //! indicate that all further calls to addSolver will be children of the current solver
  void beginChild(std::string description="");

  //! indicate the end of the current child
  void endChild();

  //! TODO update outputConnection after it has been initialize after the first transfers, or use pointer
  //! add the output connection information between two children to the current solver
  void addOutputConnection(std::shared_ptr<OutputConnection> outputConnection);

  //! set the output connector data
  template<typename FunctionSpaceType, int nComponents1, int nComponents2>
  void setOutputConnectorData(std::shared_ptr<Data::OutputConnectorData<FunctionSpaceType,nComponents1,nComponents2>> outputConnectorData);

  //! set the output connector data
  template<typename T>
  void setOutputConnectorData(std::shared_ptr<std::vector<T>> outputConnectorData);

  //! beginChild, addSolver etc. as well as writeDiagramFile will have no effect, this is needed for MultipleInstances such that not all instances appear as separate solvers but only the first one
  void disable();

  //! beginChild and addSolver and writeDiagramFile etc. will work again
  void enable();

  //! produce the resulting file, only if solverAddingEnabled_
  void writeDiagramFile(std::string filename);

  //! get the diagram as string
  std::string getDiagram();

  /** Representation of one solver in the diagram
   */
  struct solver_t
  {
    std::string name;   //< name (type) of the solver
    std::string description; //< additional string that will be included, e.g. for type of subsolver
    bool hasInternalConnectionToFirstNestedSolver;   //< if the solver has an internal connection of all output connector slots to its first subsolver. This is the case e.g. for Coupling and StrangSplitting

    /** Representation of a single output slot
     */
    struct OutputSlot
    {
      int variableNo;       //< either 0 or 1, if the slot is internally stored in variable1 or variable2 of OutputConnectorData
      std::string fieldVariableName;    //< the name of the field variable that will be written in the diagram
      std::string componentName;        //< the name of the component that will be written in the diagram
      int nComponents;      // number of components the field variable has in order to distinguish scalar field variables
    };

    std::vector<OutputSlot> outputSlots;

    std::shared_ptr<OutputConnection> outputConnection; //< pointer to the actual outputConnection object of the operator splitting

    /** connection between two output slots of two children
     */
    struct OutputConnectionRepresentation
    {
      int fromSlot;
      int toSlot;
      enum output_connection_t {ab, ba, bidirectionalCopy, bidirectionalReuse} type;    //< ab=term1 -> term2, ba=term2 -> term1, bidirectional=shared between term1 and term2
    };

    std::vector<OutputConnectionRepresentation> outputConnections;    //< connections between output slots

    std::vector<std::shared_ptr<solver_t>> children;    //< the nested solvers inside the current solver
    std::shared_ptr<solver_t> parent;                   //< pointer to the parent of the current nested solver
  };

protected:

  /** a nested class that recursively traverses the tree of nested solvers and generates the diagram
   */
  struct DiagramGenerator 
  {
    //! initialize the diagram generator, afterwards, generateDiagram can be called
    void initialize(std::shared_ptr<solver_t> solverRoot);

    //! generate the diagram
    std::string generateFinalDiagram();

  private:

    //! print the nested solver structure to result
    //! @param externalConnectionLines <fromLineNo, toLineNo, linePosition, type>
    void generateDiagramRecursion(std::stringstream &result, std::vector<std::vector<int>> &internalConnectionLinesCurrentSolver, std::vector<int> &slotLineNos, 
                                  int depth, bool isFirstChildSolver, bool isLastChildSolver);

    //! generate a diagram without the vertical connection lines, neither internal nor external lines
    std::string generateDiagramWithoutConnectionLines();
    
    //! internal variables that are used by the enclosing class to add connection lines to the diagram
    std::shared_ptr<solver_t> currentSolver_;                     //< a pointer to the current solver during the iteration
    std::vector<std::tuple<int,int,int,SolverStructureVisualizer::solver_t::OutputConnectionRepresentation::output_connection_t>> externalConnectionLines_;   //< connection lines contains the following information: <lineNoFrom, lineNoTo, lineColumn, lineType>, lineColumn is the horizontal position of the vertical data connection line
    std::vector<std::vector<int>> internalConnectionLinesAll_;    //< list of line nos that are connected by internal connections
    std::vector<int> slotLineNosAll_;                             //< all slots in the diagram with their line no
  };

  //! in the currentSolver_ fill outputConnections vector from outputConnection
  static void parseOutputConnection(std::shared_ptr<solver_t> currentSolver);

  std::shared_ptr<solver_t> solverRoot_;          //< the whole nested solver structure
  std::shared_ptr<solver_t> currentSolver_;       //< a pointer to the current solver for which call to addSolver sets the name and data

  bool enabled_;                                  //< if addSolver has an effect
  int nDisableCalls_;                             //< how often disable() has been called in sequence
};

#include "control/diagnostic_tool/solver_structure_visualizer.tpp"
