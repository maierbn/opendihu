#pragma once

#include <Python.h>  // has to be the first included header
#include <petscsys.h>
#include <list>
#include <memory>
#include <thread>
#include <map>

#ifdef HAVE_ADIOS
#include <adios2.h>
#endif

#ifdef HAVE_MEGAMOL
#include <libzmq/zmq.hpp>
#endif

#include "partition/partition_manager.h"
#include "control/python_config.h"

// forward declaration
namespace Mesh { class Manager; }
namespace Solver { class Manager; }
namespace Partition { class RankSubset; }
class SolverStructureVisualizer;

/** This class contains global variables (mesh manager and solver manager) and a python config dictionary
 *  which can be a sub-dict of the global config. Different objects of this class are used for the methods
 *  throughout this program which have a different context in terms of their config.
 */
class DihuContext
{
public:
  //! constructor, initialize context, parse command line parameters and input file
  DihuContext(int argc, char *argv[], bool doNotFinalizeMpi=false, bool settingsFromFile=true);

  //! constructor for test cases
  DihuContext(int argc, char *argv[], std::string pythonSettings, bool doNotFinalizeMpi=true);

  //! constructor with explicit pythonConfig
  DihuContext(int argc, char *argv[], bool doNotFinalizeMpi, PythonConfig pythonConfig, std::shared_ptr<Partition::RankSubset> rankSubset=nullptr);

  //! copy-constructor
  DihuContext(const DihuContext &rhs);

  //! return a context object with config originated at child node with given key
  DihuContext operator[](std::string keyString);

  //! return the top-level python config object
  PythonConfig getPythonConfig() const;

  //! return the python code that was used to create the config object
  static std::string pythonScriptText();

  //! return a text specifying the version of this opendihu program
  static std::string versionText();

  //! return a text giving meta information
  static std::string metaText();

  //! create a context object, like with the operator[] but with given config and rankSubset, if rankSubset is not given, reuse own rankSubset
  DihuContext createSubContext(PythonConfig config, std::shared_ptr<Partition::RankSubset> rankSubset=nullptr);

  //! return the mesh manager object that contains all meshes
  static std::shared_ptr<Mesh::Manager> meshManager();

  //! return the solver manager object that contains all solvers
  std::shared_ptr<Solver::Manager> solverManager() const;

  //! return the partition manager object that creates partitionings
  static std::shared_ptr<Partition::Manager> partitionManager();

  //! return the object that collects information about all nested solvers and produces a diagram, also with data connections
  static std::shared_ptr<SolverStructureVisualizer> solverStructureVisualizer();

  //! writes out the diagram file using the solverStructureVisualizer
  static void writeSolverStructureDiagram();

  //! get the own MPI rank no in the communicator of this context
  int ownRankNo();

  //! get the own MPI rank no in the world communicator
  static int ownRankNoCommWorld();

  //! number of ranks in the world communicator
  static int nRanksCommWorld();

  //! get the rank subset of this context, this may not be the same as MPI_COMM_WORLD
  std::shared_ptr<Partition::RankSubset> rankSubset() const;

#ifdef HAVE_ADIOS
  //! return the adios object
  std::shared_ptr<adios2::ADIOS> adios() const;
#endif

#ifdef HAVE_MEGAMOL
  //! get the zmq socket that can be used to send messages to MegaMol
  std::shared_ptr<zmq::socket_t> zmqSocket() const;
#endif

  //! destructor
  virtual ~DihuContext();

private:
  //! read in file and execute python script and store global variables
  bool loadPythonScriptFromFile(std::string filename);

  //! execute python script and store global variables
  void loadPythonScript(std::string text);

  //! parse the scenarioName and data under "meta" and set as global parameters
  void parseGlobalParameters();

  //! initiaize the easylogging++ framework
  void initializeLogging(int &argc, char *argv[]);

  //! start MegaMol in a new thread
  void initializeMegaMol(int argc, char *argv[]);

  //! initialize mpi context for adios
  void initializeAdios(int argc, char *argv[]);

  //! initialize python interpreter
  void initializePython(int argc, char *argv[], bool explicitConfigFileGiven);

  //! initialize the library used for network communication with MegaMol
  void initializeZMQ();

  PythonConfig pythonConfig_;                         //< the top level python config dictionary of the current context (i.e. may be a sub-dict of the global config)
  std::shared_ptr<Partition::RankSubset> rankSubset_; //< the ranks that collectively run the code where this context is valid

  static std::shared_ptr<Mesh::Manager> meshManager_; //< object that saves all meshes that are used
//  static std::shared_ptr<Solver::Manager> solverManager_; //< object that saves all solver configurations that are used
  static std::map<int, std::shared_ptr<Solver::Manager>> solverManagerForThread_;  //< object that saves all solver configurations that are used, different for each thread

  static std::shared_ptr<Partition::Manager> partitionManager_;                    //< partition manager object that creates and manages partitionings
  static std::shared_ptr<SolverStructureVisualizer> solverStructureVisualizer_;    //< object that collects information about all nested solvers and produces a diagram, also with data connections

  static int nRanksCommWorld_;                        //< number of ranks in MPI_COMM_WORLD
  static int ownRankNoCommWorld_;                     //< the own rank no in MPI_COMM_WORLD, using MPI_COMM_WORLD should be avoided in the program, instead use this global variable
  static bool initialized_;                           //< if MPI, Petsc and easyloggingPP is already initialized. This needs to be done only once in the program.
  static int nObjects_;                               //< number of objects of DihuContext, if the last object gets destroyed, call MPI_Finalize or MPI_Barrier, depending on doNotFinalizeMpi
  static std::string pythonScriptText_;               //< the text of the python config script
  static std::shared_ptr<std::thread> megamolThread_; //< thread that runs megamol
  static std::vector<char *> megamolArgv_;            //< the arguments use for the megamol instance
  static std::vector<std::string> megamolArguments_;  //< the string data of the megamol arguments
  static std::string solverStructureDiagramFile_;     //< filename of a file produced by solverStructureVisualizer_
  bool doNotFinalizeMpi_;                             //< when the last object gets destroyed, either MPI_Finalize() is called (should be used) or MPI_Barrier (only needed in testcases where MPI context needs to be used for the next test cases)

#ifdef HAVE_ADIOS
  static std::shared_ptr<adios2::ADIOS> adios_;  //< adios context option
#endif

#ifdef HAVE_MEGAMOL
  static std::shared_ptr<zmq::context_t> zmqContext_;  //< the 0mq context
  static std::shared_ptr<zmq::socket_t> zmqSocket_;  //< a socket that is connected to one megamol
#endif
};
