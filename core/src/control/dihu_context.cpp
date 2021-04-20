#include "control/dihu_context.h"

#include <Python.h>  // this has to be the first included header
#include <python_home.h>  // defines PYTHON_HOME_DIRECTORY
#include <omp.h>

#include <fstream>
#include <iostream>
#include <iomanip>
#include <memory>
#include <thread>
#include <list>
#include <petscvec.h>
#include <sys/types.h>  // getpid
#include <unistd.h>     // getpid
#include <omp.h>
#include <csignal>
#include <cstdlib>
#include <cctype>

#include "utility/python_utility.h"
//#include "output_writer/paraview/paraview.h"
#include "output_writer/python_callback/python_callback.h"
#include "output_writer/python_file/python_file.h"
#include "output_writer/exfile/exfile.h"
#include "mesh/mesh_manager/mesh_manager.h"
#include "mesh/mapping_between_meshes/manager/04_manager.h"
#include "solver/solver_manager.h"
#include "partition/partition_manager.h"
#include "control/diagnostic_tool/stimulation_logging.h"
#include "control/diagnostic_tool/solver_structure_visualizer.h"
#include "slot_connection/global_connections_by_slot_name.h"

#include "easylogging++.h"
#include "control/python_config/settings_file_name.h"
#include "utility/mpi_utility.h"
#ifdef HAVE_PAT
#include <pat_api.h>    // perftools, only available on hazel hen
#endif
#ifdef HAVE_EXTRAE
#include "extrae.h"
#endif
#ifdef HAVE_MEGAMOL
#include "Console.h"
#endif
#ifdef HAVE_CHASTE
#include "CommandLineArguments.hpp"
#include "ExecutableSupport.hpp"
#endif

bool GLOBAL_DEBUG = false;        //< use this variable to hack in debugging output that should only be visible in certain conditions. Do not commit these hacks!

// global singleton objects
std::shared_ptr<MappingBetweenMeshes::Manager>  DihuContext::mappingBetweenMeshesManager_ = nullptr;
std::shared_ptr<Mesh::Manager>                  DihuContext::meshManager_                 = nullptr;
std::shared_ptr<Solver::Manager>                DihuContext::solverManager_               = nullptr;
std::shared_ptr<Partition::Manager>             DihuContext::partitionManager_            = nullptr;
std::shared_ptr<SolverStructureVisualizer>      DihuContext::solverStructureVisualizer_   = nullptr;
std::shared_ptr<GlobalConnectionsBySlotName>    DihuContext::globalConnectionsBySlotName_ = nullptr;

// other global variables that are needed in static methods
std::string DihuContext::solverStructureDiagramFile_ = "";              //< filename of the solver structure diagram file
std::string DihuContext::pythonScriptText_ = "";                        //< the python settings text
DihuContext::logFormat_t DihuContext::logFormat_ = DihuContext::logFormat_t::logFormatCsv;    //< format of lines in the log file

// megamol variables
std::shared_ptr<std::thread> DihuContext::megamolThread_ = nullptr;
std::vector<char *> DihuContext::megamolArgv_;
std::vector<std::string> DihuContext::megamolArguments_;

#ifdef HAVE_ADIOS
std::shared_ptr<adios2::ADIOS> DihuContext::adios_ = nullptr;           //< adios context option
#endif
bool DihuContext::initialized_ = false;
int DihuContext::nObjects_ = 0;                                         //< number of objects of DihuContext, if the last object gets destroyed, call MPI_Finalize
int DihuContext::nRanksCommWorld_ = 0;                                  //< number of MPI ranks in MPI_COMM_WORLD
int DihuContext::ownRankNoCommWorld_ = 0;                               //< own MPI rank no in MPI_COMM_WORLD

void handleSignal(int signalNo)
{
  std::string signalName = strsignal(signalNo);
  Control::PerformanceMeasurement::setParameter("exit_signal",signalNo);
  Control::PerformanceMeasurement::setParameter("exit",signalName);
  Control::PerformanceMeasurement::writeLogFile();
  Control::StimulationLogging::writeLogFile();
  DihuContext::writeSolverStructureDiagram();
  MappingBetweenMeshes::Manager::writeLogFile();

  int rankNo = DihuContext::ownRankNoCommWorld();
  LOG(INFO) << "Rank " << rankNo << " received signal " << sys_siglist[signalNo]
    << " (" << signalNo << "): " << signalName;

  if (signalNo == SIGBUS)
  {
    LOG(ERROR) << "Available memory was exceeded.";
  }

  if (signalNo != SIGRTMIN)
  {
    MPI_Abort(MPI_COMM_WORLD,0);
  }

  if (signalNo == SIGSEGV)
  {
#ifndef NDEBUG
#ifndef RELWITHDEBINFO
#ifdef __GNUC__
    // source: https://stackoverflow.com/questions/77005/how-to-automatically-generate-a-stacktrace-when-my-program-crashes
    void *array[100];

    // get void*'s for all entries on the stack
    size_t size = backtrace(array, 100);

    // print stack trace
    backtrace_symbols_fd(array, size, STDERR_FILENO);
#endif
#endif
#endif
  }

  // set back to normal in case program continues execution
  Control::PerformanceMeasurement::setParameter("exit","normal");
}

// copy-constructor
DihuContext::DihuContext(const DihuContext &rhs) : pythonConfig_(rhs.pythonConfig_), rankSubset_(rhs.rankSubset_)
{
  nObjects_++;
  VLOG(1) << "DihuContext(a), nObjects = " << nObjects_;

  doNotFinalizeMpi_ = rhs.doNotFinalizeMpi_;
}

DihuContext::DihuContext(int argc, char *argv[], bool doNotFinalizeMpi, PythonConfig pythonConfig, std::shared_ptr<Partition::RankSubset> rankSubset) :
  pythonConfig_(pythonConfig), rankSubset_(rankSubset), doNotFinalizeMpi_(doNotFinalizeMpi)
{
  nObjects_++;
  VLOG(1) << "DihuContext(b), nObjects = " << nObjects_;

  // if rank subset was not given
  if (!rankSubset_)
  {
    rankSubset_ = std::make_shared<Partition::RankSubset>();   // create rankSubset with all ranks, i.e. MPI_COMM_WORLD
  }
}

DihuContext::DihuContext(int argc, char *argv[], bool doNotFinalizeMpi, bool settingsFromFile) :
  pythonConfig_(NULL), doNotFinalizeMpi_(doNotFinalizeMpi)
{
  nObjects_++;
  VLOG(1) << "DihuContext(c), nObjects = " << nObjects_;

  if (!initialized_)
  {

#ifdef HAVE_PAT
    PAT_record(PAT_STATE_OFF);
#endif

    // initialize MPI, this is necessary to be able to call PetscFinalize without MPI shutting down
    MPI_Init(&argc, &argv);

   // the following three lines output the MPI version during compilation, use for debugging
   //#define XSTR(x) STR(x)
   //#define STR(x) #x
   //#pragma message "The value of MPI_VERSION: " XSTR(MPI_VERSION)

#if MPI_VERSION >= 3
    MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN);   // MPI >= 4
#else
    MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
#endif

    // get current MPI version
    char versionString[MPI_MAX_LIBRARY_VERSION_STRING];
    int resultLength;
    MPI_Get_library_version(versionString, &resultLength);

    std::string mpiVersion(versionString, versionString+resultLength);

#ifdef HAVE_EXTRAE
    // Disable Extrae tracing at startup
    Extrae_shutdown();
#endif

    // create rankSubset with all ranks, i.e. MPI_COMM_WORLD
    rankSubset_ = std::make_shared<Partition::RankSubset>();   

    // initial global variables that will be returned by the static methods nRanksCommWorld() and ownRankNoCommWorld()
    nRanksCommWorld_ = rankSubset_->size();
    ownRankNoCommWorld_ = rankSubset_->ownRankNo();

    // initialize the logging output formats of easylogging++
    initializeLogging(argc, argv);

    // configure PETSc to abort on error
    PetscOptionsSetValue(NULL, "-on_error_abort", "");

    // initialize PETSc
    PetscInitialize(&argc, &argv, NULL, "This is an opendihu application.");

    // print header text to console
    LOG(INFO) << "This is " << versionText() << ", " << metaText();
    LOG(INFO) << mpiVersion;
    LOG(DEBUG) << "MPI version: \"" << mpiVersion << "\".";

    // warn if OpenMPI 4 is used, remove this warning if you know if the bug has been fixed (try running fibers_emg with at least 64 ranks)
    if (mpiVersion.find("Open MPI v4") != std::string::npos && metaText().find("ipvs-epyc") != std::string::npos)
    {
      LOG(WARNING) << "Using MPI 4 on ipvs-epyc, which might cause problems. \n"
        << "In 2020 we found there is a bug in at least OpenMPI 4.0.4. Everything works fine with OpenMPI 3 (e.g. version 3.1.6). \n"
        << "So, either use OpenMPI 3 or you might check if the issues have already be fixed in newer versions of OpenMPI 4 or higher. \n"
        << "If so, remove this warning message.";
    }
    // output process ID in debug
    int pid = getpid();
    LOG(DEBUG) << "PID " << pid;

    // set number of OpenMP threads to use to 1
    omp_set_num_threads(1);
    LOG(DEBUG) << "set number of threads to 1";
    
    // parallel debugging barrier
    bool enableDebuggingBarrier = false;
    PetscErrorCode ierr = PetscOptionsHasName(NULL, NULL, "-pause", (PetscBool *)&enableDebuggingBarrier); CHKERRV(ierr);

    if (enableDebuggingBarrier)
    {
      MPIUtility::gdbParallelDebuggingBarrier();
    }

    // register signal handler functions on various signals. This enforces dumping of the log file even if the program crashes.
    struct sigaction signalHandler;

    signalHandler.sa_handler = handleSignal;
    sigemptyset(&signalHandler.sa_mask);
    signalHandler.sa_flags = 0;

    sigaction(SIGINT, &signalHandler, NULL);
    //sigaction(SIGKILL, &signalHandler, NULL);
    sigaction(SIGTERM, &signalHandler, NULL);
    //sigaction(SIGABRT, &signalHandler, NULL);
    sigaction(SIGFPE, &signalHandler, NULL);
    sigaction(SIGILL, &signalHandler, NULL);
    sigaction(SIGSEGV, &signalHandler, NULL);
    sigaction(SIGXCPU, &signalHandler, NULL);
    sigaction(SIGRTMIN, &signalHandler, NULL);
    Control::PerformanceMeasurement::setParameter("exit","normal");
    Control::PerformanceMeasurement::setParameter("exit_signal","");

    // determine settings filename
    Control::settingsFileName = "settings.py";

    // check if the first command line argument is *.py, only then it is treated as config file
    bool explicitConfigFileGiven = false;
    if (argc > 1 && settingsFromFile)
    {
      std::string firstArgument = argv[1];
      if (firstArgument.rfind(".py") == firstArgument.size() - 3)
      {
        explicitConfigFileGiven = true;
        Control::settingsFileName = argv[1];
      }
      else
      {
        LOG(ERROR) << "First command line argument does not have suffix *.py, not considering it as config file!";
      }
    }

    initializePython(argc, argv, explicitConfigFileGiven);

    // load python script
    if (settingsFromFile)
    {
      if (!loadPythonScriptFromFile(Control::settingsFileName))
      {
        // if a settings file was given but could not be loaded
        if (explicitConfigFileGiven)
        {
          LOG(FATAL) << "Could not load settings file \"" << Control::settingsFileName << "\".";
        }
        else
        {
          // look for any other settings.py files to output a message with a suggestion
          std::stringstream commandSuggestions;
          int ret = system("ls ../settings*.py > a");
          if (ret == 0)
          {
            // parse contents of ls command
            std::ifstream file("a", std::ios::in|std::ios::binary);
            if (file.is_open())
            {
              while (!file.eof())
              {
                std::string line;
                std::getline(file, line);
                if (file.eof())
                  break;

                commandSuggestions << "  " << argv[0] << " " << line << std::endl;
              }
            }
            // remove temporary file `a`
            ret = system("rm a");
          }
          if (commandSuggestions.str().empty())
          {
            commandSuggestions << "  " << argv[0] << " ../settings.py";
          }

          // if no settings file was given (default file "settings.py" was tried but not successful)
          LOG(FATAL) << "No settings file was specified!" << std::endl
            << "Usually you run the executable from within a \"build_release\" or \"build_debug\" directory "
            << "and the settings file is located one directory higher. Try a command like the following:" << std::endl << std::endl
            << commandSuggestions.str();
        }
      }
    }

    // start megamol console
    LOG(DEBUG) << "initializeMegaMol";
    initializeMegaMol(argc, argv);

    // initialize chaste, if available
#ifdef HAVE_CHASTE
    int *argcChaste = new int;
    char ***argvChaste = new (char**);
    *argcChaste = argc;
    *argvChaste = argv;
    ExecutableSupport::StartupWithoutShowingCopyright(argcChaste, argvChaste);
    LOG(DEBUG) << "initialize chaste";
#endif

    initialized_ = true;
  }

  if (!rankSubset_)
    rankSubset_ = std::make_shared<Partition::RankSubset>();   // create rankSubset with all ranks, i.e. MPI_COMM_WORLD

  // if this is the first constructed DihuContext object, create global objects partition manager, mesh manager and solver manager
  if (!partitionManager_)
  {
    VLOG(2) << "create partitionManager_";
    partitionManager_ = std::make_shared<Partition::Manager>(pythonConfig_);
  }

  if (!meshManager_)
  {
    VLOG(2) << "create meshManager_";
    meshManager_ = std::make_shared<Mesh::Manager>(pythonConfig_);
    meshManager_->setPartitionManager(partitionManager_);
    mappingBetweenMeshesManager_ = std::make_shared<MappingBetweenMeshes::Manager>(pythonConfig_);
  }
  
  if (!solverManager_)
  {
    solverManager_ = std::make_shared<Solver::Manager>(pythonConfig_);
  }

  if (!solverStructureVisualizer_)
  {
    solverStructureVisualizer_ = std::make_shared<SolverStructureVisualizer>();
  }

  if (!globalConnectionsBySlotName_)
  {
    globalConnectionsBySlotName_ = std::make_shared<GlobalConnectionsBySlotName>(pythonConfig_);
  }
}

DihuContext::DihuContext(int argc, char *argv[], std::string pythonSettings, bool doNotFinalizeMpi) :
  DihuContext(argc, argv, doNotFinalizeMpi, false)
{
  nObjects_++;
  VLOG(1) << "DihuContext(d), nObjects = " << nObjects_;
  // This constructor is called when creating the context object from unit tests.

  // load python config script as given by parameter
  loadPythonScript(pythonSettings);
  if (VLOG_IS_ON(1))
  {
    PythonUtility::printDict(pythonConfig_.pyObject());
  }

  // initialize global singletons
  // they are first set to nullptr to allow break points on their destruction (happens in unit tests)
  partitionManager_ = nullptr;
  partitionManager_ = std::make_shared<Partition::Manager>(pythonConfig_);
  
  VLOG(2) << "recreate meshManager";
  meshManager_ = nullptr;
  meshManager_ = std::make_shared<Mesh::Manager>(pythonConfig_);
  meshManager_->setPartitionManager(partitionManager_);

  mappingBetweenMeshesManager_ = nullptr;
  mappingBetweenMeshesManager_ = std::make_shared<MappingBetweenMeshes::Manager>(pythonConfig_);
  
  // create solver manager
  solverManager_ = nullptr;
  solverManager_ = std::make_shared<Solver::Manager>(pythonConfig_);
  
}

PythonConfig DihuContext::getPythonConfig() const
{
  return pythonConfig_;
}

std::string DihuContext::pythonScriptText()
{
  return pythonScriptText_;
}

std::string DihuContext::versionText()
{
  std::stringstream versionTextStr;

  versionTextStr << "opendihu 1.3, built " << __DATE__; // << " " << __TIME__; // do not add time otherwise it wants to recompile this file every time
#ifdef __cplusplus
  versionTextStr << ", C++ " << __cplusplus;
#endif

#ifdef __INTEL_COMPILER
  versionTextStr << ", Intel";
#elif defined _CRAYC
  versionTextStr << ", Cray";
#elif defined __GNUC__
  versionTextStr << ", GCC";
#elif defined __PGI
  versionTextStr << ", PGI";  
#endif
#ifdef __VERSION__
  versionTextStr << " " << __VERSION__;
#elif defined __PGIC__
  versionTextStr << " " << __PGIC__;
#endif

  return versionTextStr.str();
}

std::string DihuContext::metaText()
{
  std::stringstream metaTextStr;

  // time stamp
  auto t = std::time(nullptr);
  auto tm = *std::localtime(&t);
  // metaTextStr << "current time: " << std::put_time(&tm, "%Y/%m/%d %H:%M:%S") << ", hostname: ";
  std::string tm_string = StringUtility::timeToString(&tm);
  metaTextStr << "current time: " << tm_string << ", hostname: ";

  // host name
  char hostname[MAXHOSTNAMELEN+1];
  gethostname(hostname, MAXHOSTNAMELEN+1);
  metaTextStr << std::string(hostname) << ", n ranks: " << nRanksCommWorld_;

  return metaTextStr.str();
}

int DihuContext::ownRankNoCommWorld()
{
  return ownRankNoCommWorld_;
}

int DihuContext::ownRankNo()
{
  return rankSubset_->ownRankNo();
}

int DihuContext::nRanksCommWorld()
{
  return nRanksCommWorld_;
}


std::shared_ptr<Mesh::Manager> DihuContext::meshManager()
{
  return meshManager_;
}

std::shared_ptr<MappingBetweenMeshes::Manager> DihuContext::mappingBetweenMeshesManager()
{
  return mappingBetweenMeshesManager_;
}

std::shared_ptr<Partition::Manager> DihuContext::partitionManager()
{
  return partitionManager_;
}

std::shared_ptr<Solver::Manager> DihuContext::solverManager() const
{
  return solverManager_;
}

std::shared_ptr<SolverStructureVisualizer> DihuContext::solverStructureVisualizer()
{
  return solverStructureVisualizer_;
}

std::shared_ptr<GlobalConnectionsBySlotName> DihuContext::globalConnectionsBySlotName()
{
  return globalConnectionsBySlotName_;
}

void DihuContext::writeSolverStructureDiagram()
{
  if (solverStructureVisualizer_ && solverStructureDiagramFile_ != "")
    solverStructureVisualizer_->writeDiagramFile(solverStructureDiagramFile_);
}


#ifdef HAVE_ADIOS
std::shared_ptr<adios2::ADIOS> DihuContext::adios() const
{
  return adios_;
}
#endif

#ifdef HAVE_MEGAMOL
std::shared_ptr<zmq::socket_t> DihuContext::zmqSocket() const
{
  return zmqSocket_;
}
#endif

std::shared_ptr<Partition::RankSubset> DihuContext::rankSubset() const
{
  return rankSubset_;
}

DihuContext::logFormat_t DihuContext::logFormat() {
  return logFormat_;
}

void DihuContext::setLogFormat(DihuContext::logFormat_t format) {
  logFormat_ = format;
}

DihuContext DihuContext::operator[](std::string keyString)
{
  int argc = 0;
  char **argv = NULL;
  DihuContext dihuContext(argc, argv, doNotFinalizeMpi_, PythonConfig(pythonConfig_, keyString), rankSubset_);

  return dihuContext;
}

//! create a context object, like with the operator[] but with given config
DihuContext DihuContext::createSubContext(PythonConfig config, std::shared_ptr<Partition::RankSubset> rankSubset)
{
  int argc = 0;
  char **argv = NULL;
  DihuContext dihuContext(argc, argv, doNotFinalizeMpi_, config, rankSubset);

  return dihuContext;
}

DihuContext::~DihuContext()
{
  nObjects_--;

  VLOG(1) << "~DihuContext, nObjects = " << nObjects_;
  if (nObjects_ == 0)
  {
    // write log files
    writeSolverStructureDiagram();
    Control::StimulationLogging::writeLogFile();
    Control::PerformanceMeasurement::writeLogFile();
    MappingBetweenMeshes::Manager::writeLogFile();

    // After a call to MPI_Finalize we cannot call MPI_Initialize() anymore.
    // This is only a problem when the code is tested with the GoogleTest framework, because then we want to run multiple tests in one executable.
    // In this case, do not finalize MPI, but call MPI_Barrier instead which also syncs the ranks.

    if (!doNotFinalizeMpi_)
    {
#ifdef HAVE_MEGAMOL
      LOG(DEBUG) << "wait for MegaMol to finish";

      // wait for megamol to finish
      if (megamolThread_)
        megamolThread_->join();
#endif

      // global barrier
      MPI_Barrier(MPI_COMM_WORLD);

      // finalize Petsc and MPI
      LOG(DEBUG) << "Petsc_Finalize";
      PetscErrorCode ierr = PetscFinalize(); CHKERRV(ierr);
      MPI_Finalize();
    }
  }
  // do not clear pythonConfig_ here, because it crashes
  //VLOG(4) << "PY_CLEAR(PYTHONCONFIG_)";  // note: calling VLOG in a destructor is critical and can segfault
  //Py_CLEAR(pythonConfig_);

  // do not finalize Python because otherwise tests keep crashing
  //Py_Finalize();
#if PY_MAJOR_VERSION >= 3
  //Py_Finalize();
#endif

  // do not finalize Petsc because otherwise there can't be multiple DihuContext objects for testing
  //PetscErrorCode ierr;
  //ierr = PetscFinalize(); CHKERRV(ierr);
  //MPI_Finalize();
}
