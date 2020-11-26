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
  
  std::array<PetscInt,2> globalSize_ = {5,5};
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
  std::array<PetscInt,2> beginElementGlobal_, nElementsLocal_, nRanks_
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
  
  LOG(DEBUG)
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
  // run with 2 processes
  
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
  //
  // values
  // |1.1  3.1|  |2.2  3.2|
  // |0.1  2.1|  |0.2  1.2|
  //
  // |1.1  |5.3  3.2|
  // |0.1  |2.3  1.2|

  
  // keywords: VecCreateGhost, VecGhostUpdateBegin, VecGetLocalVector
  
  MPI_Init(&argc, &argv);
  PetscInitialize(&argc, &argv, NULL, NULL);
  initializeLogging(argc, argv);
  
  PetscMPIInt ownRankNo, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &ownRankNo);
  MPI_Comm_size(PETSC_COMM_WORLD,&size);
  
  LOG(DEBUG) << " ownRankNo: " << ownRankNo << ", size: " << size;

  PetscInt nNodesLocal = 0;      // non-ghosts
  PetscInt nNodesGlobal = 0;
  PetscInt nGhosts = 0;
  
  PetscErrorCode ierr;
  Vec globalVector;
  Vec localVector;
  ISLocalToGlobalMapping localToGlobalMapping;
  
  if (ownRankNo == 0)
  {
    nNodesLocal = 2;   // non-ghosts
    nNodesGlobal = 6;
    nGhosts = 2;
    std::array<PetscInt,2> ghostDofGlobalNos({2,4});
 
    ierr = VecCreateGhost(MPI_COMM_WORLD, nNodesLocal, nNodesGlobal, nGhosts, ghostDofGlobalNos.data(), &globalVector); CHKERRQ(ierr);
    VecZeroEntries(globalVector);
    
    ierr = VecGetLocalToGlobalMapping(globalVector, &localToGlobalMapping); CHKERRQ(ierr);
  
    /*VecGetLocalVector(globalVector, localVector);
    VecRestoreLocalVector(globalVector, localVector);*/
    
    VecGhostUpdateBegin(globalVector, INSERT_VALUES, SCATTER_FORWARD);
    VecGhostUpdateEnd(globalVector, INSERT_VALUES, SCATTER_FORWARD);
    
    VecGhostGetLocalForm(globalVector, &localVector);
    
    std::array<PetscInt,4> indices({0,1,2,3});
    std::array<double,4> values({0.1,1.1,2.1,3.1});
      
    // here VecSetValues is used, not VecSetValuesLocal, because it is a plain vector
    VecSetValues(localVector, 4, indices.data(), values.data(), INSERT_VALUES);
    
    VecGhostRestoreLocalForm(globalVector, &localVector);
    
    VecGhostUpdateBegin(globalVector, ADD_VALUES, SCATTER_REVERSE);
    VecGhostUpdateEnd(globalVector, ADD_VALUES, SCATTER_REVERSE);
    /*
    VecGetLocalVectorRead(globalVector, localVector);
    VecRestoreLocalVectorRead(globalVector, localVector);*/
    //VecGhostGetLocalForm(Vec g,Vec *l)//
    
  }
  else if (ownRankNo == 1)
  {
    nNodesLocal = 4;
    nNodesGlobal = 6;
    nGhosts = 0;
    std::array<PetscInt,0> ghostDofGlobalNos;
 
    ierr = VecCreateGhost(MPI_COMM_WORLD, nNodesLocal, nNodesGlobal, nGhosts, ghostDofGlobalNos.data(), &globalVector); CHKERRQ(ierr);
    VecZeroEntries(globalVector);
    
    ierr = VecGetLocalToGlobalMapping(globalVector, &localToGlobalMapping); CHKERRQ(ierr);
  
    VecGhostUpdateBegin(globalVector, INSERT_VALUES, SCATTER_FORWARD);
    VecGhostUpdateEnd(globalVector, INSERT_VALUES, SCATTER_FORWARD);
    
    VecGhostGetLocalForm(globalVector, &localVector);
    
    std::array<PetscInt,4> indices({0,1,2,3});
    std::array<double,4> values({0.2,1.2,2.2,3.2});
      
    VecSetValues(localVector, 4, indices.data(), values.data(), INSERT_VALUES);
    
    VecGhostRestoreLocalForm(globalVector, &localVector);
    
    VecGhostUpdateBegin(globalVector, ADD_VALUES, SCATTER_REVERSE);
    VecGhostUpdateEnd(globalVector, ADD_VALUES, SCATTER_REVERSE);
    
    /*VecGetLocalVector(globalVector, localVector);
    VecRestoreLocalVector(globalVector, localVector);
    
    VecGetLocalVectorRead(globalVector, localVector);
    VecRestoreLocalVectorRead(globalVector, localVector);*/
  }
  
  LOG(DEBUG) << "localToGlobalMapping: " << localToGlobalMapping;

  int nEntriesGlobal;
  VecGetSize(globalVector, &nEntriesGlobal);  

  LOG(DEBUG) << "global vector has " << nEntriesGlobal << " entries";
  
  std::vector<int> indices(nEntriesGlobal);
  std::iota(indices.begin(), indices.end(), 0);
  std::vector<double> values(nEntriesGlobal);
    
  VecGetValues(globalVector, nEntriesGlobal, indices.data(), values.data());
  LOG(DEBUG) << "global values: " << values;
  
  int nEntriesLocal;
  VecGetSize(localVector, &nEntriesLocal);
  LOG(DEBUG) << "local vector has " << nEntriesLocal << " entries ("
    << nNodesLocal << " local, " << (nGhosts) << " ghosts)";
  
  std::vector<int> indices2(nEntriesLocal);
  std::iota(indices2.begin(), indices2.end(), 0);
  std::vector<double> values2(nEntriesLocal);
    
  VecGetValues(localVector, nEntriesLocal, indices2.data(), values2.data());
  
  LOG(DEBUG) << "local values: " << values2;
  VecView(globalVector, PETSC_VIEWER_STDOUT_WORLD);

  MPI_Barrier(MPI_COMM_WORLD);
  LOG(INFO) << "communicate ghost values";

  // now communicate ghost values from rank 1 to rank 0
  VecGhostUpdateBegin(globalVector, INSERT_VALUES, SCATTER_FORWARD);
  VecGhostUpdateEnd(globalVector, INSERT_VALUES, SCATTER_FORWARD);

  VecGetValues(globalVector, nEntriesGlobal, indices.data(), values.data());
  LOG(INFO) << "global values: " << values;

  VecGetValues(localVector, nEntriesLocal, indices2.data(), values2.data());
  LOG(INFO) << "local values: " << values2;

  if (ownRankNo == 0)
  {
    // global values
    assert(fabs(values[0] - 0.1) < 1e-12);
    assert(fabs(values[1] - 1.1) < 1e-12);
    assert(fabs(values[2] - 2.3) < 1e-12);
    assert(fabs(values[3] - 5.3) < 1e-12);

    // local values
    assert(fabs(values2[0] - 0.1) < 1e-12);
    assert(fabs(values2[1] - 1.1) < 1e-12);
    assert(fabs(values2[2] - 2.3) < 1e-12);
    assert(fabs(values2[3] - 5.3) < 1e-12);
  }
  else
  {
    // global values
    assert(fabs(values[2] - 2.3) < 1e-12);
    assert(fabs(values[3] - 1.2) < 1e-12);
    assert(fabs(values[4] - 5.3) < 1e-12);
    assert(fabs(values[5] - 3.2) < 1e-12);

    // local values
    assert(fabs(values2[0] - 2.3) < 1e-12);
    assert(fabs(values2[1] - 1.2) < 1e-12);
    assert(fabs(values2[2] - 5.3) < 1e-12);
    assert(fabs(values2[3] - 3.2) < 1e-12);
  }


  MPI_Barrier(MPI_COMM_WORLD);
  // --------- matrix ------------
  LOG(DEBUG) << " --------- matrix ---------";
  Mat globalMatrix;
  Mat localMatrix;
  int nComponents = 1;
  
  ierr = MatCreate(MPI_COMM_WORLD, &globalMatrix); CHKERRQ(ierr);
  //ierr = MatSetType(globalMatrix, MATMAIJ); CHKERRQ(ierr);
  ierr = MatSetFromOptions(globalMatrix); CHKERRQ(ierr);
  ierr = MatSetSizes(globalMatrix, nComponents*nNodesLocal, nComponents*nNodesLocal, nComponents*nNodesGlobal, nComponents*nNodesGlobal); CHKERRQ(ierr);
  ierr = MatSetUp(globalMatrix);  // either MatSetUp or MatSetPreallocation has to be used
  ierr = MatSetLocalToGlobalMapping(globalMatrix, localToGlobalMapping, localToGlobalMapping); CHKERRQ(ierr);
    
  
  int nRows, nColumns, nRowsLocal, nColumnsLocal;
  ierr = MatGetSize(globalMatrix, &nRows, &nColumns); CHKERRQ(ierr);
  ierr = MatGetLocalSize(globalMatrix, &nRowsLocal, &nColumnsLocal); CHKERRQ(ierr);
  LOG(DEBUG) << "matrix global: " << nRows << "x" << nColumns << ", local: " << nRowsLocal << "x" << nColumnsLocal;

  std::vector<int> localDofNos(nNodesLocal+nGhosts);
  std::iota(localDofNos.begin(), localDofNos.end(), 0);
  
  LOG(DEBUG) << "matrix created";
  
  IS indexSet;
  ierr = ISCreateGeneral(PETSC_COMM_WORLD, localDofNos.size(), localDofNos.data(), PETSC_COPY_VALUES, &indexSet); CHKERRQ(ierr);
  ierr = MatGetLocalSubMatrix(globalMatrix, indexSet, indexSet, &localMatrix); CHKERRQ(ierr);
  // see Figure 10 on p.72 of PETSc manual
  
  LOG(DEBUG) << "got local submatrix for indexSet " << localDofNos;
  
  // set values
  std::vector<int> rowIndices;
  std::vector<int> columnIndices;
  
  //nGhosts = 0;
  rowIndices.resize(nNodesLocal);
  std::iota(rowIndices.begin(), rowIndices.end(), 0);
  
  columnIndices.resize(nNodesLocal);
  std::iota(columnIndices.begin(), columnIndices.end(), 0);
  
  values.resize(std::pow((nNodesLocal+nGhosts),2));
  for (int i = 0; i < std::pow((nNodesLocal+nGhosts),2); i++)
  {
    values[i] = 1.1*(ownRankNo+1);
  }
  LOG(DEBUG) << "INSERT_VALUES: rowIndices: " << rowIndices << ", columnIndices: " << columnIndices << ", values " << values;
  
  ierr = MatSetValuesLocal(localMatrix, rowIndices.size(), rowIndices.data(), columnIndices.size(), columnIndices.data(), values.data(), INSERT_VALUES); CHKERRQ(ierr);
  //ierr = MatSetValuesLocal(globalMatrix, rowIndices.size(), rowIndices.data(), columnIndices.size(), columnIndices.data(), values.data(), ADD_VALUES); CHKERRQ(ierr);
  
  ierr = MatRestoreLocalSubMatrix(globalMatrix, indexSet, indexSet, &localMatrix); CHKERRQ(ierr);
  
  ierr = MatAssemblyBegin(globalMatrix, MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(globalMatrix, MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
  
  
  ierr = MatGetLocalSubMatrix(globalMatrix, indexSet, indexSet, &localMatrix); CHKERRQ(ierr);
  
  LOG(DEBUG) << "ADD_VALUES: rowIndices: " << rowIndices << ", columnIndices: " << columnIndices << ", values " << values;
  
  rowIndices.resize(nNodesLocal+nGhosts);
  std::iota(rowIndices.begin(), rowIndices.end(), 0);
  
  columnIndices.resize(nNodesLocal+nGhosts);
  std::iota(columnIndices.begin(), columnIndices.end(), 0);
  
  ierr = MatSetValuesLocal(localMatrix, rowIndices.size(), rowIndices.data(), columnIndices.size(), columnIndices.data(), values.data(), ADD_VALUES); CHKERRQ(ierr);
  
  
  // test set values
  rowIndices.resize(1);
  rowIndices[0] = 0;
  columnIndices.resize(1);
  columnIndices[0] = 0;
  values[0] = 10.0;
  //ierr = MatSetValuesLocal(localMatrix, rowIndices.size(), rowIndices.data(), columnIndices.size(), columnIndices.data(), values.data(), INSERT_VALUES); CHKERRQ(ierr);
  
  ierr = MatRestoreLocalSubMatrix(globalMatrix, indexSet, indexSet, &localMatrix); CHKERRQ(ierr);
  
  ierr = MatAssemblyBegin(globalMatrix, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(globalMatrix, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  
//  ierr = MatAssemblyBegin(globalMatrix, MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
//  ierr = MatAssemblyEnd(globalMatrix, MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
  
  // show entries
  ierr = MatView(globalMatrix, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  //exit(0);
  // again get submatrix
  ierr = MatGetLocalSubMatrix(globalMatrix, indexSet, indexSet, &localMatrix); CHKERRQ(ierr);
  
  rowIndices.resize(1);
  rowIndices[0] = 1;
  const int diag = 77;
  ierr = MatZeroRowsColumnsLocal(globalMatrix, rowIndices.size(), rowIndices.data(), diag, NULL, NULL); CHKERRQ(ierr);
  
  //ierr = MatRestoreLocalSubMatrix(globalMatrix, indexSet, indexSet, &localMatrix); CHKERRQ(ierr);
  
  
  ierr = MatAssemblyBegin(globalMatrix, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(globalMatrix, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  
  MPI_Barrier(MPI_COMM_WORLD);
  // get all values 
  //LOG(DEBUG) << " display globalMatrix";
  //LOG(DEBUG) << " globalMatrix: " << PetscUtility::getStringMatrix(globalMatrix);
  
  // show entries
  ierr = MatView(globalMatrix, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  PetscFinalize();
  MPI_Finalize();
  
  return EXIT_SUCCESS;
}
