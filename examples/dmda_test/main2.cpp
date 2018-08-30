#include <iostream>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <array>

#include "easylogging++.h"
#include "semt/Semt.h"

#include "opendihu.h"
#include <petscdmda.h>
/*
void testDmda()
{
  PetscErrorCode ierr;
  
  std::array<int,2> globalSize_ = {5,5};
  int nDofsPerElement = 1;
  int ghostLayerWidth = 1;
  
  DM da;
  
  // create 2d decomposition
  ierr = DMDACreate2d(MPI_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DMDA_STENCIL_BOX,
                      globalSize_[0], globalSize_[1], PETSC_DECIDE, PETSC_DECIDE,
                      nDofsPerElement, ghostLayerWidth, NULL, NULL, &da); CHKERRV(ierr);
                      
  // get global coordinates of local partition
  PetscInt x, y, m, n;
  ierr = DMDAGetCorners(da, &x, &y, NULL, &m, &n, NULL); CHKERRV(ierr);
  std::array<int,2> beginElementGlobal_, nElementsLocal_, nRanks_
  std::array<std::vector<int>,2> localSizesOnRanks_;
  beginElementGlobal_[0] = (global_no_t)x;
  beginElementGlobal_[1] = (global_no_t)y;
  nElementsLocal_[0] = (element_no_t)m;
  nElementsLocal_[1] = (element_no_t)n;
  
  // get number of ranks in each coordinate direction
  ierr = DMDAGetInfo(da, NULL, NULL, NULL, NULL, &nRanks_[0], &nRanks_[1], NULL, NULL, NULL, NULL, NULL, NULL, NULL); CHKERRV(ierr);
  
  // get local sizes on the ranks
  const PetscInt *lxData;
  const PetscInt *lyData;
  ierr = DMDAGetOwnershipRanges(da, &lxData, &lyData, NULL);
  localSizesOnRanks_[0].assign(lxData, lxData + nRanks_[0]);
  localSizesOnRanks_[1].assign(lyData, lyData + nRanks_[1]);
  
  LOG(INFO)  
    << "globalSize_: " << globalSize_
    << "beginElementGlobal_: " << beginElementGlobal_ << ", "
    << "nElementsLocal_: " << nElementsLocal_ << ", "
    << "nRanks_: " << nRanks_ << ", "
    << "localSizesOnRanks_: " << localSizesOnRanks_;
  
}
*/

void initializeLogging(int argc, char *argv[])
{
  START_EASYLOGGINGPP(argc, argv);
/*
  std::ifstream file("logging.conf");
  if (!file.is_open())
  {
    // if file does not exist, create it
    std::ofstream out("logging.conf");
    if (!out.is_open())
    {
      LOG(ERROR) << "Could not open logging file for output";
    }
    out << R"(
* GLOBAL:
   FORMAT               =  "INFO : %msg"
   FILENAME             =  "/tmp/logs/my.log"
   ENABLED              =  true
   TO_FILE              =  true
   TO_STANDARD_OUTPUT   =  true
   SUBSECOND_PRECISION  =  1
   PERFORMANCE_TRACKING =  false
   MAX_LOG_FILE_SIZE    =  2097152 ## 2MB - Comment starts with two hashes (##)
   LOG_FLUSH_THRESHOLD  =  100 ## Flush after every 100 logs
* DEBUG:
   FORMAT               = "DEBUG: %msg"
* WARNING:
   FORMAT               = "WARN : %loc %func: Warning: %msg"
* ERROR:
   FORMAT               = "ERROR: %loc %func: Error: %msg"
* FATAL:
   FORMAT               = "FATAL: %loc %func: Fatal error: %msg"
    )";
  }
  file.close();

  el::Configurations conf("logging.conf");
*/

// color codes: https://github.com/shiena/ansicolor/blob/master/README.md
#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_LIGHT_GRAY    "\x1b[90m"
#define ANSI_COLOR_LIGHT_WHITE    "\x1b[97m"
#define ANSI_COLOR_RESET   "\x1b[0m"

  std::string separator(80, '_');
  el::Configurations conf;
  conf.setToDefault();

  int nRanks;
  int rankNo;
  MPI_Comm_size(MPI_COMM_WORLD, &nRanks);
  MPI_Comm_rank(MPI_COMM_WORLD, &rankNo);
  
  std::string prefix;
  if (nRanks > 1)
  {
    std::stringstream s;
    s << rankNo << "/" << nRanks << " ";
    prefix = s.str();
  }
  
  conf.setGlobally(el::ConfigurationType::Format, prefix+"INFO : %msg");
  conf.setGlobally(el::ConfigurationType::Filename, "/tmp/logs/opendihu.log");
  conf.setGlobally(el::ConfigurationType::Enabled, "true");
  conf.setGlobally(el::ConfigurationType::ToFile, "true");
  conf.setGlobally(el::ConfigurationType::ToStandardOutput, "true");

  // set format of outputs
  conf.set(el::Level::Debug, el::ConfigurationType::Format, prefix+"DEBUG: %msg");
  conf.set(el::Level::Trace, el::ConfigurationType::Format, prefix+"TRACE: %msg");
  conf.set(el::Level::Verbose, el::ConfigurationType::Format, ANSI_COLOR_LIGHT_WHITE "" + prefix+"VERB%vlevel: %msg" ANSI_COLOR_RESET);
  conf.set(el::Level::Warning, el::ConfigurationType::Format,
           prefix+"WARN : %loc %func: \n" ANSI_COLOR_YELLOW "Warning: " ANSI_COLOR_RESET "%msg");

  conf.set(el::Level::Error, el::ConfigurationType::Format,
           prefix+"ERROR: %loc %func: \n" ANSI_COLOR_RED "Error: %msg" ANSI_COLOR_RESET);

  conf.set(el::Level::Fatal, el::ConfigurationType::Format,
           std::string(ANSI_COLOR_MAGENTA)+prefix+"FATAL: %loc %func: \n"+separator
           +"\nFatal error: %msg\n"+separator+ANSI_COLOR_RESET+"\n");

  //el::Loggers::addFlag(el::LoggingFlag::HierarchicalLogging);

//#ifdef NDEBUG      // if release
//  conf.set(el::Level::Debug, el::ConfigurationType::Enabled, "false");
//  std::cout<<"DISABLE Debug"<<std::endl;
//#endif

  // reconfigure all loggers
  el::Loggers::reconfigureAllLoggers(conf);
  el::Loggers::removeFlag(el::LoggingFlag::AllowVerboseIfModuleNotSpecified);
}

int main(int argc, char *argv[])
{
  
  // test case
  
  // global
  // |3 |4 5|
  // |0 |1 2|
  //  p0 p1
  
  // petsc
  // |1 |4 5|
  // |0 |2 3|
  //  p0 p1
  
  // local
  // |1 3*| |2 3|
  // |0 2*| |0 1|
  //  p0     p1   *=ghost
  
  
  // keywords: VecCreateGhost, VecGhostUpdateBegin, VecGetLocalVector
  
  MPI_Init(&argc, &argv);
  PetscInitialize(&argc, &argv, NULL, NULL);
  initializeLogging(argc, argv);
  
  PetscMPIInt ownRankNo, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &ownRankNo);
  MPI_Comm_size(PETSC_COMM_WORLD,&size);
  
  LOG(INFO) << " ownRankNo: " << ownRankNo << ", size: " << size;

  int nNodesLocal = 2;      // non-ghosts
  int nNodesGlobal = 4;
  int nGhosts = 0;
  
  PetscErrorCode ierr;
  // --------- matrix ------------
  LOG(INFO) << " --------- matrix ---------";
  Mat globalMatrix;
  //Mat localMatrix;
  int nComponents = 1;
  
  ierr = MatCreate(MPI_COMM_WORLD, &globalMatrix); CHKERRQ(ierr);
  ierr = MatSetType(globalMatrix, MATMAIJ); CHKERRQ(ierr);
  ierr = MatSetFromOptions(globalMatrix); CHKERRQ(ierr);
  //ierr = MatSetSizes(globalMatrix, nComponents*nNodesLocal, nComponents*nNodesLocal, nComponents*nNodesGlobal, nComponents*nNodesGlobal); CHKERRQ(ierr);
  ierr = MatSetSizes(globalMatrix, PETSC_DECIDE, PETSC_DECIDE, nComponents*nNodesGlobal, nComponents*nNodesGlobal); CHKERRQ(ierr);

  int diagonalNonZeros = nNodesGlobal;
  int offdiagonalNonZeros = nNodesGlobal;
  ierr = MatMPIAIJSetPreallocation(globalMatrix, diagonalNonZeros, NULL, offdiagonalNonZeros, NULL); CHKERRQ(ierr);
  ierr = MatSeqAIJSetPreallocation(globalMatrix, diagonalNonZeros, NULL); CHKERRQ(ierr);
    
  
  int nRows, nColumns, nRowsLocal, nColumnsLocal;
  ierr = MatGetSize(globalMatrix, &nRows, &nColumns); CHKERRQ(ierr);
  ierr = MatGetLocalSize(globalMatrix, &nRowsLocal, &nColumnsLocal); CHKERRQ(ierr);
  LOG(INFO) << "matrix global: " << nRows << "x" << nColumns << ", local: " << nRowsLocal << "x" << nColumnsLocal;

  
  //ierr = MatSetLocalToGlobalMapping(globalMatrix, localToGlobalMapping, localToGlobalMapping); CHKERRQ(ierr);
  
  //MatSetValues(globalMatrix, 1);
  /*
  std::vector<int> localDofNos(nNodesLocal);
  std::iota(localDofNos.begin(), localDofNos.end(), 0);
  
  LOG(INFO) << "matrix created";
  
  IS indexSet;
  ierr = ISCreateGeneral(PETSC_COMM_SELF, localDofNos.size(), localDofNos.data(), PETSC_COPY_VALUES, &indexSet); CHKERRQ(ierr);
  */
  //ierr = MatGetLocalSubMatrix(globalMatrix, indexSet, indexSet, &localMatrix); CHKERRQ(ierr);
  // see Figure 10 on p.72 of PETSc manual
  
  //LOG(INFO) << "got local submatrix for indexSet " << localDofNos;
  
  // set values
  std::vector<int> rowIndices;
  std::vector<int> columnIndices;
  
  rowIndices.resize(nNodesLocal+nGhosts);
  std::iota(rowIndices.begin(), rowIndices.end(), 0);
  
  columnIndices.resize(nNodesLocal+nGhosts);
  std::iota(columnIndices.begin(), columnIndices.end(), 0);
  
  std::vector<double> values;
  values.resize(std::pow((nNodesLocal),2));
  for (int i = 0; i < std::pow((nNodesLocal),2); i++)
  {
    values[i] = 1.1*(ownRankNo+1);
  }
  
  LOG(INFO) << "rowIndices: " << rowIndices << ", columnIndices: " << columnIndices << ", set values " << values;
  
  //ierr = MatSetValuesLocal(localMatrix, rowIndices.size(), rowIndices.data(), columnIndices.size(), columnIndices.data(), values.data(), ADD_VALUES); CHKERRQ(ierr);
  //ierr = MatSetValuesLocal(globalMatrix, rowIndices.size(), rowIndices.data(), columnIndices.size(), columnIndices.data(), values.data(), ADD_VALUES); CHKERRQ(ierr);
  ierr = MatSetValues(globalMatrix, rowIndices.size(), rowIndices.data(), columnIndices.size(), columnIndices.data(), values.data(), ADD_VALUES); CHKERRQ(ierr);
  
  //LOG(INFO) << "restore local submatrix";
  //ierr = MatRestoreLocalSubMatrix(globalMatrix, indexSet, indexSet, &localMatrix); CHKERRQ(ierr);
  
  
  ierr = MatAssemblyBegin(globalMatrix, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(globalMatrix, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  
  //ierr = MatAssemblyBegin(globalMatrix, MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
  //ierr = MatAssemblyEnd(globalMatrix, MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
  
  
  // again get submatrix
  //ierr = MatGetLocalSubMatrix(globalMatrix, indexSet, indexSet, &localMatrix); CHKERRQ(ierr);
    
  //MPI_Barrier(MPI_COMM_WORLD);
  //LOG(INFO) << " display localMatrix";
  //ierr = MatAssemblyBegin(localMatrix, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  //ierr = MatAssemblyEnd(localMatrix, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  //LOG(INFO) << " localMatrix: " << PetscUtility::getStringMatrix(localMatrix);
  
  //ierr = MatRestoreLocalSubMatrix(globalMatrix, indexSet, indexSet, &localMatrix); CHKERRQ(ierr);
  
  
  MPI_Barrier(MPI_COMM_WORLD);
  // get all values 
  //LOG(INFO) << " display globalMatrix";
  //LOG(INFO) << " globalMatrix: " << PetscUtility::getStringMatrix(globalMatrix);
  
  // show entries
  ierr = MatView(globalMatrix, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  
  PetscFinalize();
  MPI_Finalize();
  
  return EXIT_SUCCESS;
}