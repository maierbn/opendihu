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

#include "utility/python_utility.h"
#include "output_writer/paraview/paraview.h"
#include "output_writer/python_callback/python_callback.h"
#include "output_writer/python_file/python_file.h"
#include "output_writer/exfile/exfile.h"
#include "mesh/mesh_manager.h"
#include "solver/solver_manager.h"
#include "partition/partition_manager.h"

#include "easylogging++.h"
#include "control/settings_file_name.h"
#include "utility/mpi_utility.h"
#ifdef HAVE_PAT
#include <pat_api.h>    // perftools, only available on hazel hen
#endif
#ifdef HAVE_MEGAMOL
#include "Console.h"
#endif

//INITIALIZE_EASYLOGGINGPP

std::shared_ptr<Mesh::Manager> DihuContext::meshManager_ = nullptr;
//std::shared_ptr<Solver::Manager> DihuContext::solverManager_ = nullptr;
std::map<int, std::shared_ptr<Solver::Manager>> DihuContext::solverManagerForThread_;
std::shared_ptr<Partition::Manager> DihuContext::partitionManager_ = nullptr;
std::shared_ptr<std::thread> DihuContext::megamolThread_ = nullptr;
std::vector<char *> DihuContext::megamolArgv_;
std::vector<std::string> DihuContext::megamolArguments_;

#ifdef HAVE_ADIOS
std::shared_ptr<adios2::ADIOS> DihuContext::adios_ = nullptr;  ///< adios context option
std::shared_ptr<adios2::IO> DihuContext::io_ = nullptr;        ///< IO object of adios
#endif
bool DihuContext::initialized_ = false;
int DihuContext::nObjects_ = 0;   ///< number of objects of DihuContext, if the last object gets destroyed, call MPI_Finalize
int DihuContext::nRanksCommWorld_ = 0;   ///< number of objects of DihuContext, if the last object gets destroyed, call MPI_Finalize

// copy-constructor
DihuContext::DihuContext(const DihuContext &rhs) : pythonConfig_(rhs.pythonConfig_)
{
  nObjects_++;
  VLOG(1) << "DihuContext(a), nObjects = " << nObjects_;

  doNotFinalizeMpi_ = rhs.doNotFinalizeMpi_;
}


DihuContext::DihuContext(int argc, char *argv[], bool doNotFinalizeMpi, PythonConfig pythonConfig) :
  pythonConfig_(pythonConfig), doNotFinalizeMpi_(doNotFinalizeMpi)
{
  nObjects_++;
  VLOG(1) << "DihuContext(b), nObjects = " << nObjects_;
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

    // get global number of MPI ranks
    MPIUtility::handleReturnValue (MPI_Comm_size(MPI_COMM_WORLD, &nRanksCommWorld_));

    // load configuration from file if it exits
    initializeLogging(argc, argv);

    // configure PETSc to abort on errorm
    PetscOptionsSetValue(NULL, "-on_error_abort", "");

    // initialize PETSc
    PetscInitialize(&argc, &argv, NULL, "This is an opendihu application.");

    // set number of threads to use to 1
    omp_set_num_threads(1);
    LOG(DEBUG) << "set number of threads to 1";

    // output process ID
    int pid = getpid();
    LOG(DEBUG) << "PID " << pid;

    // parallel debugging barrier
    bool enableDebuggingBarrier = false;
    PetscErrorCode ierr = PetscOptionsHasName(NULL, NULL, "-pause", (PetscBool *)&enableDebuggingBarrier); CHKERRV(ierr);

    if (enableDebuggingBarrier)
    {
      MPIUtility::gdbParallelDebuggingBarrier();
    }

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
      loadPythonScriptFromFile(Control::settingsFileName);
    }

    // start megamol console
    LOG(DEBUG) << "initializeMegaMol";
    initializeMegaMol(argc, argv);

    initialized_ = true;
  }
  
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
  }
  
  if (solverManagerForThread_.empty())
  {
    VLOG(2) << "create solverManagerForThread_";
    // create solver manager for thread 0
    solverManagerForThread_[0] = std::make_shared<Solver::Manager>(pythonConfig_);
    solverManagerForThread_[1] = std::make_shared<Solver::Manager>(pythonConfig_);
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

  partitionManager_ = nullptr;
  partitionManager_ = std::make_shared<Partition::Manager>(pythonConfig_);
  
  VLOG(2) << "recreate meshManager";
  meshManager_ = nullptr;
  meshManager_ = std::make_shared<Mesh::Manager>(pythonConfig_);
  meshManager_->setPartitionManager(partitionManager_);
  
  // create solver manager for thread 0
  solverManagerForThread_.clear();
  solverManagerForThread_[0] = std::make_shared<Solver::Manager>(pythonConfig_);
  
}

PythonConfig DihuContext::getPythonConfig() const
{
  return pythonConfig_;
}

int DihuContext::ownRankNo()
{
  int rankNo;
  MPIUtility::handleReturnValue (MPI_Comm_rank(MPI_COMM_WORLD, &rankNo));
  return rankNo;
}

std::shared_ptr<Mesh::Manager> DihuContext::meshManager()
{
  return meshManager_;
}

std::shared_ptr<Partition::Manager> DihuContext::partitionManager()
{
  return partitionManager_;
}

std::shared_ptr<Solver::Manager> DihuContext::solverManager() const
{
  // get number of omp threads
  //int nThreads = omp_get_num_threads();
  int threadId = omp_get_thread_num();
  
  if (solverManagerForThread_.find(threadId) == solverManagerForThread_.end())
  {
    LOG(DEBUG) << "create solver manager for thread " << threadId;
    // create solver manager
    solverManagerForThread_[threadId] = std::make_shared<Solver::Manager>(pythonConfig_);
    
    LOG(DEBUG) << "(done)";
  }
  else 
  {
    LOG(DEBUG) << "solver manager for thread " << threadId << " exists";
  }
  
  return solverManagerForThread_[threadId];
}

#ifdef HAVE_ADIOS
std::shared_ptr<adios2::IO> DihuContext::adiosIo()
{
  return io_;
}
#endif

DihuContext DihuContext::operator[](std::string keyString)
{
  int argc = 0;
  char **argv = NULL;
  DihuContext dihuContext(argc, argv, doNotFinalizeMpi_, PythonConfig(pythonConfig_, keyString));

  return dihuContext;
}

//! create a context object, like with the operator[] but with given config
DihuContext DihuContext::createSubContext(PythonConfig config)
{
  int argc = 0;
  char **argv = NULL;
  DihuContext dihuContext(argc, argv, doNotFinalizeMpi_, config);

  return dihuContext;
}

DihuContext::~DihuContext()
{
  nObjects_--;

  VLOG(1) << "~DihuContext, nObjects = " << nObjects_;
  if (nObjects_ == 0)
  {
    // write log file
    Control::PerformanceMeasurement::writeLogFile();

    // After a call to MPI_Finalize we cannot call MPI_Initialize() anymore.
    // This is only a problem when the code is tested with the GoogleTest framework, because then we want to run multiple tests in one executable.
    // In this case, do not finalize MPI, but call MPI_Barrier instead which also syncs the ranks.

    if (doNotFinalizeMpi_)
    {
      //LOG(DEBUG) << "MPI_Barrier";
      //MPI_Barrier(MPI_COMM_WORLD);
    }
    else
    {
#ifdef HAVE_MEGAMOL
      LOG(DEBUG) << "wait for MegaMol to finish";

      // wait for megamol to finish
      megamolThread_->join();
#endif

      LOG(DEBUG) << "MPI_Finalize";
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
