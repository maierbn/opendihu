// parallel-in-time (XBraid) Implicit Euler Solver ()

#include "specialized_solver/parallel_in_time/MultiDomain/PinT_MD.h"

#include <omp.h>
#include <sstream>

#include <braid.h>
#include "braid_test.h"
#include "specialized_solver/parallel_in_time/MultiDomain/PinT_MD_Braid.h"
#include "specialized_solver/parallel_in_time/MultiDomain/PinT_lib_MD.h"
#include "specialized_solver/parallel_in_time/MultiDomain/PinT_fun_MD.h"

#include <petscdm.h>
#include <petscdmda.h>
#include <petscts.h>
#include <petscdraw.h>
#include <petscvec.h>

namespace ParallelInTime
{

template<class NestedSolverMD>
PinTMD<NestedSolverMD>::
PinTMD(DihuContext context) :
  Runnable(),
  context_(context["PinTMD"]), initialized_(false)
{
  // get python settings object from context
  this->specificSettings_ = this->context_.getPythonConfig();

  // initialize output writers
  this->outputWriterManager_.initialize(this->context_, this->specificSettings_);

  // parse options
  int myOption = this->specificSettings_.getOptionInt("myOption", 1, PythonUtility::Positive);

  LOG(DEBUG) << "myOption: " << myOption;
}


template<class NestedSolverMD>
void PinTMD<NestedSolverMD>::
initialize()
{
  // make sure that we initialize only once, in the next call, initialized_ is true
  if (initialized_)
    return;
  /*
  std::vector<PyObject *> MultiDomainConfigs;

  PyObject *MultiDomainConfig = this->specificSettings_.template getOptionListBegin<PyObject *>("TimeSteppingScheme");

  // loop over other entries of list
  for (;
    !this->specificSettings_.getOptionListEnd("TimeSteppingScheme");
    this->specificSettings_.template getOptionListNext<PyObject *>("TimeSteppingScheme", MultiDomainConfig))
  {
    MultiDomainConfigs.push_back(MultiDomainConfig);
  }
  */
  //pid_t pid = getpid();
  //printf("pid: %d", pid);

  //std::this_thread::sleep_for (std::chrono::seconds(30));

  // parse all settings in TimeSteppingScheme
  PyObject *MultiDomainConfig = this->specificSettings_.template getOptionListBegin<PyObject *>("TimeSteppingScheme");
  for (;
  !this->specificSettings_.getOptionListEnd("TimeSteppingScheme");
  this->specificSettings_.template getOptionListNext<PyObject *>("TimeSteppingScheme", MultiDomainConfig)){}

  // split MPI communicator, create communicators with rank numbering in x domain
  // MPI_Comm communicatorTotal = MPI_COMM_WORLD;

  int nRanksInSpace = this->specificSettings_.getOptionInt("nRanksInSpace", 1, PythonUtility::Positive);
  braid_SplitCommworld(&communicatorTotal_, nRanksInSpace, &communicatorX_, &communicatorT_);
  
  //PinT_initialize();

  LOG(DEBUG) << "before creation of rank subset X";

  // create rankSubset and assign to partitionManager, to be used by all further created meshes and solvers
  rankSubsetX_ = std::make_shared<Partition::RankSubset>(communicatorX_);
  DihuContext::partitionManager()->setRankSubsetForNextCreatedPartitioning(rankSubsetX_);

  LOG(DEBUG) << "nextRankSubset: " << *DihuContext::partitionManager()->nextRankSubset();

  // int size;
  // MPI_Comm_size(communicatorX, &size);
  //
  // LOG(DEBUG) << "rankSubsetX: " << size;
  // LOG(DEBUG) << "rankSubsetX: " << *rankSubsetX_;

  //int test;
  //test = 0;
  // loop over parsed config objects of implicit eulers
  /*
  for (int i = 0; i < MultiDomainConfigs.size(); i++)
  {
    LOG(DEBUG) << "i = " << i;

    PyObject *MultiDomainConfig = MultiDomainConfigs[i];
  */

    LOG(DEBUG) << "create sub context";

    //LOG(DEBUG) << "MultiDomainConfig: " << MultiDomainConfig;
    DihuContext MultiDomainContext = context_.createSubContext(MultiDomainConfig, rankSubsetX_);

    LOG(DEBUG) << "nextRankSubset: " << *DihuContext::partitionManager()->nextRankSubset();

    //LOG(DEBUG) << "MultiDomainContext: " << MultiDomainContext.getPythonConfig();

    //if (test == 0) {
    //  if (i==0){
    //    braid_SplitCommworld(&communicatorTotal_, 1, &communicatorX, &communicatorT);
    //
    //    // create rankSubset and assign to partitionManager, to be used by all further created meshes and solvers
    //    rankSubsetX_ = std::make_shared<Partition::RankSubset>(communicatorX);
    //    DihuContext::partitionManager()->setRankSubsetForNextCreatedPartitioning(rankSubsetX_);
    //    test=1;
    //  }
    //  else {
    //    // create rank subset
    //    braid_SplitCommworld(&communicatorTotal_, nRanksInSpace, &communicatorX, &communicatorT);
    //
    //    // create rankSubset and assign to partitionManager, to be used by all further created meshes and solvers
    //    rankSubsetX_ = std::make_shared<Partition::RankSubset>(communicatorX);
    //    DihuContext::partitionManager()->setRankSubsetForNextCreatedPartitioning(rankSubsetX_);
    //  }
    //}
    LOG(DEBUG) << "create new MultidomainSolver";
    // create rank subset
    //std::shared_ptr<Partition::RankSubset> nextRankSubset = std::make_shared<Partition::RankSubset>(communicatorX_);
    ///DihuContext::partitionManager()->setRankSubsetForNextCreatedPartitioning(nextRankSubset);
  
    MultiDomainSolvers_.push_back(
      std::make_shared<NestedSolverMD>(MultiDomainContext)
    );

    LOG(DEBUG) << "nextRankSubset: " << *DihuContext::partitionManager()->nextRankSubset();

    LOG(DEBUG) << "puskh bvack multidomain context";

    data_.push_back(
      std::make_shared<Data>(MultiDomainContext)
    );

    //LOG(DEBUG) << "initialize multidomain for i = " << i;

    LOG(DEBUG) << "initialize MultidomainSolvers";

    MultiDomainSolvers_.back()->initialize();
  

  //MultiDomainSolvers_[0](std::make_shared<NestedSolverMD>(MultiDomainContext));
  //MultiDomainSolvers_[0](std::make_shared<Data>(MultiDomainContext));
  //MultiDomainSolvers_[0]->initialize();

    LOG(DEBUG) << "retrieve functionSpace";
    // In order to initialize the data object and actuall create all variables, we first need to assign a function space to the data object.
    // A function space object of type FunctionSpace<MeshType,BasisFunctionType> (see "function_space/function_space.h")
    // is an object that stores the mesh (e.g., all nodes and elements) as well as the basis function (e.g. linear Lagrange basis functions).
    // The NestedSolverMD solver already created a function space that we should use. We already have a typedef "FunctionSpace" that is the class of NestedSolverMD's function space type.
    std::shared_ptr<FunctionSpace> functionSpace = MultiDomainSolvers_.back()->data().functionSpace();
    //std::shared_ptr<FunctionSpace> functionSpace = MultiDomainSolvers_[0]->data().functionSpace();

    // Pass the function space to the data object. data_ stores field variables.
    // It needs to know the number of nodes and degrees of freedom (dof) from the function space in order to create the vectors with the right size.
    data_.back()->setFunctionSpace(functionSpace);
    //data_[0]->setFunctionSpace(functionSpace);

    LOG(DEBUG) << "initialize data";

    //LOG(DEBUG) << "initialize data for i = " << i;
    // now call initialize, data will then create all variables (Petsc Vec's)
    data_.back()->initialize();
    //data_[0]->initialize();


    // it is also possible to pass some field variables from the data of the NestedSolverMD to own data object
    // data_.back()->setSolutionVariable(MultiDomainSolvers_.back()->data().solution());
  //};


  LOG(DEBUG) << "n: " << MultiDomainSolvers_.size();

  // initialize PinT
  PinT_initialize();
  // here is the space to initialize anything else that is needed for your solver
  // for example, initialize Braid here (not like this)
  initialized_ = true;
}

template<class NestedSolverMD>
void PinTMD<NestedSolverMD>::
run()
{


  // initialize
  initialize();

  int          rank;
  double       loglevels;
  PetscInt     nspace        =  this->nspace_+1;
  nspace = this->MultiDomainSolvers_[0]->nSolutionValuesLocal();


  // Define XBraid parameters

  // int       max_levels    = 3;
  int       nrelax        = 3;
  int       skip          = 0;
  double    tol           = 1.0e-07;
  int       cfactor       = 2;
  int       max_iter      = 100;
  int       min_coarse    = 3;
  int       fmg           = 0;
  int       scoarsen      = 0;
  int       res           = 0;
  int       wrapper_tests = 0;
  // int       print_level   = 2;
  int       access_level  = 1;
  int       use_sequential= 0;

  // communicatorTotal   = MPI_COMM_WORLD;
  MPI_Comm_rank(communicatorTotal_, &rank);

  /* The first step before running simulations, is always to verify the wrapper tests */
  if(wrapper_tests)
  {
     /* Create spatial communicator for wrapper-tests */
     /*braid_SplitCommworld(&comm, 1, &comm_x, &comm_t);*/

     braid_TestAll(app_, communicatorX_, stdout, 0.0, (tstop_-tstart_)/ntime_,
                   2*(tstop_-tstart_)/ntime_, my_Init_MD, my_Free_MD, my_Clone_MD,
                   my_Sum_MD, my_SpatialNorm_MD, my_BufSize_MD, my_BufPack_MD,
                   my_BufUnpack_MD, my_Coarsen_MD, my_Interp_MD, my_Residual_MD, my_Step_MD);
  }
  else
  {
     /* Scale tol by domain */
     tol = tol/( sqrt((tstop_ - tstart_)/(ntime_-1))*sqrt((xstop_ - xstart_)/(nspace-1)) );

     /* Set Braid options */
     braid_SetPrintLevel( core_, print_level_);
     braid_SetAccessLevel( core_, access_level);
     braid_SetMaxLevels(core_, max_levels_);
     braid_SetMinCoarse( core_, min_coarse );
     braid_SetSkip(core_, skip);
     braid_SetNRelax(core_, -1, nrelax);
     braid_SetAbsTol(core_, tol);
     braid_SetCFactor(core_, -1, cfactor);
     braid_SetMaxIter(core_, max_iter);
     braid_SetSeqSoln(core_, use_sequential);
     if (fmg)
     {
        braid_SetFMG(core_);
     }
     if (res)
     {
        braid_SetResidual(core_, my_Residual_MD);
     }
     loglevels = log2(nspace - 1.0);
     if ( scoarsen && ( fabs(loglevels - round(loglevels)) > 1e-10 ))
     {
        if(rank == 0)
        {
           fprintf(stderr, "\nWarning!\nFor spatial coarsening, spatial grids must be a "
                   "power of 2 + 1, \ni.e., nspace = 2^k + 1.  Your spatial grid is of size"
                   " %d.  Spatial \ncoarsening is therefore being ignored.\n\n", nspace);
        }
     }
     else if (scoarsen)
     {
        braid_SetSpatialCoarsen(core_, my_Coarsen_MD);
        braid_SetSpatialRefine(core_,  my_Interp_MD);
     }

     LOG(DEBUG) << "........................";

     braid_Drive(core_);

     /* Print accumulated info on space-time grids visited during the simulation */
     if( (print_level_ > 0) && (rank == 0))
     {
        print_sc_info(app_->sc_info, max_levels_);
     }
  }
  /* Clean up */
  braid_Destroy(core_);
  free( app_->sc_info);
  //free( app_->g);
  free( app_ );

  // do something else
  //executeMyHelperMethod();

  // write current output values using the output writers
  //this->outputWriterManager_.writeOutput(*this->data_.back());
}

template<class NestedSolverMD>
void PinTMD<NestedSolverMD>::
reset()
{
  initialized_ = false;
  // "uninitialize" everything
}

template<class NestedSolverMD>
void PinTMD<NestedSolverMD>::
executeMyHelperMethod()
{/*
  // this is the template for an own private method

  // Get the solution field variable from the data object. The Petsc object can be retrieved with "valuesGlobal()".
  Vec solution = data_.solution()->valuesGlobal();

  // As an example we invert all the solution values.
  // Because "Vec"'s are actually pointers, this effects the actual data, no copy-back is needed.
  // Note the error handling with Petsc functions. Always use "ierr = PetscFunction(); CHKERRV(ierr);"
  PetscErrorCode ierr;
  ierr = VecScale(solution, -1); CHKERRV(ierr);

  // get a field variable from data object:
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace,1>> fieldVariableB = data_.fieldVariableB();

  // copy the solution from the NestedSolverMD_ to fieldVariableB.
  ierr = VecCopy(solution, fieldVariableB->valuesGlobal()); CHKERRV(ierr);

  // add 1.0 to every entry in fieldVariableA, this already updates fieldVariableA in data because it is a pointer. There is no need to copy the values back.
  ierr = VecShift(fieldVariableB->valuesGlobal(), 1.0); CHKERRV(ierr);

  // // For example you can also get the stiffness matrix of the NestedSolverMD_
  // std::shared_ptr<PartitionedPetscMat<FunctionSpace>> stiffnessMatrix = NestedSolverMD_.data().stiffnessMatrix();
  // Mat m = stiffnessMatrix->valuesGlobal();
  //
  // // e.g. add identity to m
  // ierr = MatShift(m, 1.0); CHKERRV(ierr);*/
}

template<class NestedSolverMD>
void PinTMD<NestedSolverMD>::
PinT_initialize()
{
  // initialize time stepping values
  // tstart_ = 0.0;
  // tstop_ = 1.0;
  // ntime_ = 10;
  // nspace_=8;
  // PetscReal *initialGuess_=[2,2,4,5,2,2];
  if (specificSettings_.hasKey("tstart"))
    tstart_ = specificSettings_.getOptionDouble("tstart", 0.0);
  if (specificSettings_.hasKey("tstop"))
    tstop_ = specificSettings_.getOptionDouble("tstop", 1.0, PythonUtility::Positive);
  if (specificSettings_.hasKey("ntime"))
    ntime_ = specificSettings_.getOptionDouble("ntime", 1.0, PythonUtility::Positive);
  if (specificSettings_.hasKey("nspace"))
    nspace_ = specificSettings_.getOptionDouble("nspace", 1.0, PythonUtility::Positive);

  int       nspace        = this->nspace_+1;
  nspace = this->MultiDomainSolvers_[0]->nSolutionValuesLocal();

  app_ = (my_App *) malloc(sizeof(my_App));
  //(app_->g)             = (double*) malloc( nspace*sizeof(double) );
  (app_->comm)          = communicatorTotal_;
  (app_->tstart)        = tstart_;
  (app_->tstop)         = tstop_;
  (app_->ntime)         = ntime_;
  (app_->xstart)        = xstart_;
  (app_->xstop)         = xstop_;
  (app_->nspace)        = nspace;
  (app_->print_level)   = print_level_;
  (app_->MultiDomainSolvers)        = &this->MultiDomainSolvers_;

  /* Initialize storage for sc_info, for tracking space-time grids visited during the simulation */
  app_->sc_info = (double*) malloc( 2*max_levels_*sizeof(double) );
  for(int i = 0; i < 2*max_levels_; i++) {
     app_->sc_info[i] = -1.0;
  }

  /* Initialize Braid */
  braid_Init(communicatorTotal_, communicatorT_, tstart_, tstop_, ntime_, app_,
         my_Step_MD, my_Init_MD, my_Clone_MD, my_Free_MD, my_Sum_MD, my_SpatialNorm_MD,
         my_Access_MD, my_BufSize_MD, my_BufPack_MD, my_BufUnpack_MD, &core_);

}

template<class NestedSolverMD>
typename PinTMD<NestedSolverMD>::Data &PinTMD<NestedSolverMD>::
data()
{
  // get a reference to the data object
  return *data_.back();

  // The NestedSolverMD_ object also has a data object, we could also directly use this and avoid having an own data object:
  //  return NestedSolverMD_.data();
}

//! get the data that will be transferred in the operator splitting to the other term of the splitting
//! the transfer is done by the output_connector_data_transfer class
template<class NestedSolverMD>
std::shared_ptr<typename PinTMD<NestedSolverMD>::OutputConnectorDataType> PinTMD<NestedSolverMD>::
getOutputConnectorData()
{
  //! This is relevant only, if this solver is part of a splitting or coupling scheme. Then this method returns the values/variables that will be
  // transferred to the other solvers. We can just reuse the values of the NestedSolverMD_.
  return MultiDomainSolvers_.back()->getOutputConnectorData();
}

// template<class NestedSolverMD>
// void PinTMD<NestedSolverMD>::
// setTimeSpan(double tstart, double tstop)
// {
  
//   (app_->tstart)        = tstart;
//   (app_->tstop)         = tstop;

//   //if (timeStepWidth_ > endTime_-startTime_)
//   //{
//   //  LOG(DEBUG) << "time span [" << startTime << "," << endTime << "], reduce timeStepWidth from " << timeStepWidth_ << " to " << endTime_-startTime_;
//   //  timeStepWidth_ = endTime_-startTime_;
//   //}

//   /* Initialize Braid */
//   braid_Init(MPI_COMM_WORLD, communicatorTotal_, tstart_, tstop_, ntime_, app_,
//          my_Step_MD, my_Init_MD, my_Clone, my_Free, my_Sum, my_SpatialNorm,
//          my_Access, my_BufSize, my_BufPack, my_BufUnpack, &core_);

// }

// template<class NestedSolverMD>
// void PinTMD<NestedSolverMD>::
// advanceTimeSpan()
// {
//   run();
// }


} //namespace ParallelInTime
