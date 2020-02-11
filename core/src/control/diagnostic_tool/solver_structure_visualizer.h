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
  void addSolver(std::string name);

  //! indicate that all further calls to addSolver will be children of the current solver
  void beginChild(std::string description="");

  //! indicate the end of the current child
  void endChild();

  //! add the output connection information between two children to the current solver
  void addOutputConnection(OutputConnection &outputConnection);

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

  /** Representation of one solver, will be a box in the SVG file
   */
  struct solver_t
  {
    std::string name;   //< name (type) of the solver
    std::string description; //< additional string that will be included, e.g. for type of subsolver

    /** Representation of a single output slot
     */
    struct OutputSlot
    {
      std::string fieldVariableName;
      std::string componentName;
      int nComponents;      // number of components the field variable has in order to distinguish scalar field variables
    };

    std::vector<OutputSlot> outputSlots;

    /** connection between two output slots of two children
     */
    struct OutputConnection
    {
      int fromSlot;
      int toSlot;
      enum output_connection_t {ab, ba, bidirectionalCopy, bidirectionalReuse} type;    //< ab=term1 -> term2, ba=term2 -> term1, bidirectional=shared between term1 and term2
    };

    std::vector<OutputConnection> outputConnections; //< connections between output slots

    std::vector<std::shared_ptr<solver_t>> children;   //< the nested solvers inside the current solver
    std::shared_ptr<solver_t> parent;                  //< pointer to the parent of the current nested solver
  };

protected:

  //! print the nested solver structure to result
  //! @param connectionLines <fromLineNo, toLineNo, linePosition, type>
  void generateDiagram(std::stringstream &result, std::vector<std::tuple<int,int,int,SolverStructureVisualizer::solver_t::OutputConnection::output_connection_t>> &connectionLines,
                       std::vector<int> &slotLineNos, int depth=0, bool isFirstChildSolver=false, bool isLastChildSolver=false);

  std::shared_ptr<solver_t> solverRoot_;          //< the whole nested solver structure
  std::shared_ptr<solver_t> currentSolver_;       //< a pointer to the current solver for which call to addSolver sets the name and data

  bool enabled_;                                  //< if addSolver has an effect
  int nDisableCalls_;                             //< how often disable() has been called in sequence
};

#include "control/diagnostic_tool/solver_structure_visualizer.tpp"
