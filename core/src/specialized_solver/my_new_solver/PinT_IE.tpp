#include "specialized_solver/my_new_solver/PinT_IE.h"

#include <omp.h>
#include <sstream>

#include "braid.h"

#include <petscdm.h>
#include <petscdmda.h>
#include <petscts.h>
#include <petscdraw.h>
#include <petscvec.h>
//#include "/mnt/c/Users/mariu/OneDrive/Dokumente/Masterarbeit/Opendihu/opendihu/dependencies/xbraid/src/xbraid-2.3.0/braid/braid.hpp"
// #include "PinT_IE_Braid.cpp"
#include "PinT_IE_Braid.c"

template<class NestedSolver>
PinT<NestedSolver>::
PinT(DihuContext context) :
  Runnable(),
  context_(context["PinT"]), initialized_(false)
{
  // get python settings object from context
  this->specificSettings_ = this->context_.getPythonConfig();

  // initialize output writers
  this->outputWriterManager_.initialize(this->context_, this->specificSettings_);

  // parse options
  int myOption = this->specificSettings_.getOptionInt("myOption", 1, PythonUtility::Positive);

  LOG(DEBUG) << "myOption: " << myOption;
}

template<class NestedSolver>
void PinT<NestedSolver>::
initialize()
{
  // make sure that we initialize only once, in the next call, initialized_ is true
  if (initialized_)
    return;

  std::vector<PyObject *> implicitEulerConfigs;

  PyObject *implicitEulerConfig = this->specificSettings_.getOptionListBegin<PyObject *>("TimeSteppingScheme");

  // loop over other entries of list
  for (;
    !this->specificSettings_.getOptionListEnd("TimeSteppingScheme");
    this->specificSettings_.getOptionListNext<PyObject *>("TimeSteppingScheme", implicitEulerConfig))
  {
    implicitEulerConfigs.push_back(implicitEulerConfig);
  }

  for (int i = 0; i < implicitEulerConfigs.size(); i++)
  {
    LOG(DEBUG) << "i = " << i;

    PyObject *implicitEulerConfig = implicitEulerConfigs[i];

    LOG(DEBUG) << "implicitEulerConfig: " << implicitEulerConfig;
    DihuContext implicitEulerContext = context_.createSubContext(implicitEulerConfig);

    LOG(DEBUG) << "implicitEulerContext: " << implicitEulerContext.getPythonConfig();


    // create rank subset
    //std::shared_ptr<RankSubset> nextRankSubet = std::make_shared<RankSubset>();

    //template<typename Iter>
    //RankSubset(Iter ranksBegin, Iter ranksEnd, std::shared_ptr<RankSubset> parentRankSubset = nullptr);

    //DihuContext::partitionManager()->setRankSubsetForNextCreatedPartitioning(nextRankSubset);


    implicitEulerSolvers_.push_back(
      std::make_shared<NestedSolver>(implicitEulerContext)
    );

    data_.push_back(
      std::make_shared<Data>(implicitEulerContext)
    );

    LOG(DEBUG) << "initialize implicit EUler for i = " << i;

    implicitEulerSolvers_.back()->initialize();

    // In order to initialize the data object and actuall create all variables, we first need to assign a function space to the data object.
    // A function space object of type FunctionSpace<MeshType,BasisFunctionType> (see "function_space/function_space.h")
    // is an object that stores the mesh (e.g., all nodes and elements) as well as the basis function (e.g. linear Lagrange basis functions).
    // The NestedSolver solver already created a function space that we should use. We already have a typedef "FunctionSpace" that is the class of NestedSolver's function space type.
    std::shared_ptr<FunctionSpace> functionSpace = implicitEulerSolvers_.back()->data().functionSpace();

    // Pass the function space to the data object. data_ stores field variables.
    // It needs to know the number of nodes and degrees of freedom (dof) from the function space in order to create the vectors with the right size.
    data_.back()->setFunctionSpace(functionSpace);


    LOG(DEBUG) << "initialize data for i = " << i;
    // now call initialize, data will then create all variables (Petsc Vec's)
    data_.back()->initialize();

    // it is also possible to pass some field variables from the data of the NestedSolver to own data object
    data_.back()->setSolutionVariable(implicitEulerSolvers_.back()->data().solution());
  };


  LOG(DEBUG) << "n: " << implicitEulerSolvers_.size();

  // here is the space to initialize anything else that is needed for your solver
  // for example, initialize Braid here (not like this)
  initialized_ = true;
}

template<class NestedSolver>
void PinT<NestedSolver>::
run()
{
  // initialize everything
  initialize();
  PinT_initialize();

  braid_Core    core;
  my_App       *app;
  // MPI_Comm      comm, comm_x, comm_t;
  MPI_Comm comm;
  // int           i, rank, arg_index;
  int           i, rank;
  double        loglevels;

  PetscReal    tstart        =  this->tstart_;
  PetscReal    tstop         =  this->tstop_;
  PetscInt     ntime         =  this->ntime_;
  PetscReal    xstart        =  0.0;
  PetscReal    xstop         =  4;
  PetscInt     nspace        =  33;

  /* Define XBraid parameters
   * See -help message for descriptions */
  int       max_levels    = 5;
  int       nrelax        = 1;
  int       skip          = 0;
  double    tol           = 1.0e-07;
  int       cfactor       = 2;
  int       max_iter      = 30;
  int       min_coarse    = 3;
  int       fmg           = 0;
  int       scoarsen      = 1;
  int       res           = 0;
  int       wrapper_tests = 0;
  int       print_level   = 3;
  int       access_level  = 1;
  int       use_sequential= 0;

  comm   = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &rank);

  app = (my_App *) malloc(sizeof(my_App));
  (app->g)             = (double*) malloc( nspace*sizeof(double) );
  (app->comm)          = comm;
  (app->tstart)        = tstart;
  (app->tstop)         = tstop;
  (app->ntime)         = ntime;
  (app->xstart)        = xstart;
  (app->xstop)         = xstop;
  (app->nspace)        = nspace;
  (app->print_level)   = print_level;
  (app->implicitEulerSolvers)        = &this->implicitEulerSolvers_;

  /* Initialize storage for sc_info, for tracking space-time grids visited during the simulation */
  app->sc_info = (double*) malloc( 2*max_levels*sizeof(double) );
  for( i = 0; i < 2*max_levels; i++) {
     app->sc_info[i] = -1.0;
  }

  /* Initialize Braid */
  braid_Init(MPI_COMM_WORLD, comm, tstart, tstop, ntime, app,
         my_Step, my_Init, my_Clone, my_Free, my_Sum, my_SpatialNorm,
         my_Access, my_BufSize, my_BufPack, my_BufUnpack, &core);

  /* The first step before running simulations, is always to verify the wrapper tests */
  if(wrapper_tests)
  {
     /* Create spatial communicator for wrapper-tests */
     /*braid_SplitCommworld(&comm, 1, &comm_x, &comm_t);*/

     /*braid_TestAll(app, comm_x, stdout, 0.0, (tstop-tstart)/ntime,
                   2*(tstop-tstart)/ntime, my_Init, my_Free, my_Clone,
                   my_Sum, my_SpatialNorm, my_BufSize, my_BufPack,
                   my_BufUnpack, my_Coarsen, my_Interp, my_Residual, my_Step); */
  }
  else
  {
     /* Scale tol by domain */
     tol = tol/( sqrt((tstop - tstart)/(ntime-1))*sqrt((xstop - xstart)/(nspace-1)) );

     /* Set Braid options */
     braid_SetPrintLevel( core, print_level);
     braid_SetAccessLevel( core, access_level);
     braid_SetMaxLevels(core, max_levels);
     braid_SetMinCoarse( core, min_coarse );
     braid_SetSkip(core, skip);
     braid_SetNRelax(core, -1, nrelax);
     braid_SetAbsTol(core, tol);
     braid_SetCFactor(core, -1, cfactor);
     braid_SetMaxIter(core, max_iter);
     braid_SetSeqSoln(core, use_sequential);
     if (fmg)
     {
        braid_SetFMG(core);
     }
     if (res)
     {
        braid_SetResidual(core, my_Residual);
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
        braid_SetSpatialCoarsen(core, my_Coarsen);
        braid_SetSpatialRefine(core,  my_Interp);
     }

     braid_Drive(core);

     /* Print accumulated info on space-time grids visited during the simulation */
     if( (print_level > 0) && (rank == 0))
     {
        print_sc_info(app->sc_info, max_levels);
     }
  }
  /* Clean up */
  braid_Destroy(core);
  free( app->sc_info);
  free( app->g);
  free( app );

  /* Finalize MPI */
  // MPI_Finalize();

  // hässlicher Versuch das Ergebnis zu vergleichen, noch nicht fertig

  // call the nested solver
  // Vec test1 = implicitEulerSolvers_[5]->data().solution()->valuesGlobal();
  // Vec U; //create vector V with the given values of u
  // PetscReal blub[33];
  // blub[0]=2;
  // blub[1]=2;
  // blub[2]=4;
  // blub[3]=5;
  // blub[4]=2;
  // blub[5]=2;
  // blub[6]=2;
  // blub[7]=2;
  // blub[8]=2;
  // blub[9]=2;
  // blub[10]=2;
  // blub[11]=2;
  // blub[12]=2;
  // blub[13]=2;
  // blub[14]=2;
  // blub[15]=2;
  // blub[16]=2;
  // blub[17]=2;
  // blub[18]=2;
  // blub[19]=2;
  // blub[20]=2;
  // blub[21]=2;
  // blub[22]=2;
  // blub[23]=2;
  // blub[24]=2;
  // blub[25]=2;
  // blub[26]=2;
  // blub[27]=2;
  // blub[28]=2;
  // blub[29]=2;
  // blub[30]=2;
  // blub[31]=2;
  // blub[32]=2;
  //
  // PetscInt numbers[]={1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32};
  // VecDuplicate(test1, &U);
  // VecSetValues(U,33,numbers,blub,INSERT_VALUES);
  // VecAssemblyBegin(U);
  // VecAssemblyEnd(U);
  // VecView(test1, 	PETSC_VIEWER_STDOUT_SELF);
  // // VecView(implicitEulerSolvers_[5]->data().solution()->valuesGlobal(), 	PETSC_VIEWER_STDOUT_SELF);
  // implicitEulerSolvers_[5]->data().solution()->valuesGlobal()=U;
  // // VecView(implicitEulerSolvers_[5]->data().solution()->valuesGlobal(), 	PETSC_VIEWER_STDOUT_SELF);
  // implicitEulerSolvers_[5]->setTimeSpan(0.0, 5.0);
  // implicitEulerSolvers_[5]->run();
  // Vec test2 = implicitEulerSolvers_[5]->data().solution()->valuesGlobal();
  // VecView(test2, 	PETSC_VIEWER_STDOUT_SELF);
  // // VecView(implicitEulerSolvers_[5]->data().solution()->valuesGlobal(), 	PETSC_VIEWER_STDOUT_SELF);
  // VecAXPY(test2, -1, test1);
  // VecView(test2, 	PETSC_VIEWER_STDOUT_SELF);
  
  // do something else
  //executeMyHelperMethod();

  // write current output values using the output writers
  this->outputWriterManager_.writeOutput(*this->data_.back());
}

template<class NestedSolver>
void PinT<NestedSolver>::
reset()
{
  initialized_ = false;
  // "uninitialize" everything
}

template<class NestedSolver>
void PinT<NestedSolver>::
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

  // copy the solution from the nestedSolver_ to fieldVariableB.
  ierr = VecCopy(solution, fieldVariableB->valuesGlobal()); CHKERRV(ierr);

  // add 1.0 to every entry in fieldVariableA, this already updates fieldVariableA in data because it is a pointer. There is no need to copy the values back.
  ierr = VecShift(fieldVariableB->valuesGlobal(), 1.0); CHKERRV(ierr);

  // // For example you can also get the stiffness matrix of the nestedSolver_
  // std::shared_ptr<PartitionedPetscMat<FunctionSpace>> stiffnessMatrix = nestedSolver_.data().stiffnessMatrix();
  // Mat m = stiffnessMatrix->valuesGlobal();
  //
  // // e.g. add identity to m
  // ierr = MatShift(m, 1.0); CHKERRV(ierr);*/
}

template<class NestedSolver>
void PinT<NestedSolver>::
PinT_initialize()
{
  // initialize time stepping values
  tstart_ = 0.0;
  tstop_ = 1.0;
  ntime_ = 10;
  // PetscReal *initialGuess_=[2,2,4,5,2,2];
  if (specificSettings_.hasKey("tstart"))
    tstart_ = specificSettings_.getOptionDouble("tstart", 0.0);
  if (specificSettings_.hasKey("tstop"))
    tstop_ = specificSettings_.getOptionDouble("tstop", 1.0, PythonUtility::Positive);
  if (specificSettings_.hasKey("ntime"))
    ntime_ = specificSettings_.getOptionDouble("ntime", 1.0, PythonUtility::Positive);
  // if (specificSettings_.hasKey("Initial Guess"))
  //   specificSettings_.getOptionVector("Initial Guess", initialGuess_);
}

template<class NestedSolver>
typename PinT<NestedSolver>::Data &PinT<NestedSolver>::
data()
{
  // get a reference to the data object
  return data_.back();

  // The nestedSolver_ object also has a data object, we could also directly use this and avoid having an own data object:
  //  return nestedSolver_.data();
}

//! get the data that will be transferred in the operator splitting to the other term of the splitting
//! the transfer is done by the output_connector_data_transfer class
template<class NestedSolver>
std::shared_ptr<typename PinT<NestedSolver>::OutputConnectorDataType> PinT<NestedSolver>::
getOutputConnectorData()
{
  //! This is relevant only, if this solver is part of a splitting or coupling scheme. Then this method returns the values/variables that will be
  // transferred to the other solvers. We can just reuse the values of the nestedSolver_.
  return implicitEulerSolvers_.back()->getOutputConnectorData();
}
