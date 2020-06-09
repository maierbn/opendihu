// specific functions for XBraid with Implicit Euler
// defines my_Step (integration step)

#include "specialized_solver/parallel_in_time/MultiDomain/PinT_MD_Braid.h"

#include "specialized_solver/parallel_in_time/MultiDomain/PinT_fun_MD.h"
#include "specialized_solver/parallel_in_time/MultiDomain/PinT_lib_MD.h"
#include "data_management/specialized_solver/PinT_MD.h"

int my_Step_MD(braid_App        app,
            braid_Vector     ustop,
            braid_Vector     fstop,
            braid_Vector     u,
            braid_StepStatus status)
{

  /*
  // use the following methods of the wrapper, instead of VecSetValues and VecGetValues

  // get the number of local solution values
  int nSolutionValuesLocal = MultiDomainSolver.nSolutionValuesLocal();

  // allocate a vector to hold all solution values (do this only once, for efficiency!)
  std::vector<double> values(nSolutionValuesLocal);

  // get the current solution values (at the end of my_Step)
  MultiDomainSolver.getSolution(values.data());
  LOG(INFO) << "number of values: " << nSolutionValuesLocal << ", values: " << values;

  // set the current solution values (at the beginning of my_Step)
  MultiDomainSolver.setSolution(values.data());

  */

   PetscReal tstart;             /* current time */
   PetscReal tstop;              /* evolve to this time*/
   PetscInt level, i, solver;
   PetscReal deltaX, deltaT;
   // PetscReal * help;

   LOG(DEBUG) << "Braid Step!";
   LOG(DEBUG) << "---------------";
   // get level and start and stop time from braid status
   braid_StepStatusGetLevel(status, &level);
   braid_StepStatusGetTstartTstop(status, &tstart, &tstop);

   // determine, which solver is needed (depends on size)
   solver=log2(u->size - 1);
   solver=0;
   assert(solver < (*app->MultiDomainSolvers).size());

   LOG(DEBUG) << "solver" << solver;
   LOG(DEBUG) << "(*app->MultiDomainSolvers).size()" << (*app->MultiDomainSolvers).size();

   /* XBraid forcing */

   if(fstop != NULL)
   {
      for(i = 0; i < u->size; i++)
      {
         u->values[i] = u->values[i] + fstop->values[i];
      }
   }

   // get alias variables
   std::shared_ptr<typename _braid_App_struct::NestedSolverMD> MultiDomainSolver = (*app->MultiDomainSolvers)[solver];
   //std::shared_ptr<typename Data::PinTMD<typename _braid_App_struct::NestedSolverMD::FunctionSpace>::ScalarFieldVariableType> solution = MultiDomainSolver->data().solution();
   // Data data;
   // auto data = MultiDomainSolver->data();
   // assert(u->size == solution->nDofsGlobal());

   LOG(DEBUG) << "u->size" << u->size;
   LOG(DEBUG) << "MultiDomainSolver->nSolutionValuesLocal()" << MultiDomainSolver->nSolutionValuesLocal();
   //int istart, iend, iterator;
   //VecGetOwnershipRange(solution->valuesGlobal(),&istart,&iend);
   assert(u->size == MultiDomainSolver->nSolutionValuesLocal());

   // Set the initial guess for the solver
   // PetscErrorCode ierr;
   // ierr = VecSetValues(solution->valuesGlobal(), u->size, solution->functionSpace()->meshPartition()->dofNosLocal().data(), u->values, INSERT_VALUES); CHKERRQ(ierr);
   // VecAssemblyBegin(solution->valuesGlobal());
   // VecAssemblyEnd(solution->valuesGlobal());
   MultiDomainSolver->setSolution(u->values);

   // set time span for the solver, which is calculated by braid status
   MultiDomainSolver->setTimeSpan(tstart, tstop);
   MultiDomainSolver->setNumberTimeSteps(10);
   MultiDomainSolver->setSystemMatrix((tstop-tstart)/1);

   // Debug Options
   LOG(DEBUG) << "--------------------------------------------------------------";
   LOG(DEBUG) << "level: " << level << ", solver: " << solver << ", size: " << u->size << ", t: [" << tstart << "," << tstop << "], before implicit euler:" << u->values;
   
   // PetscRealView(u->size, u->values, 0);
   // VecView(MultiDomainSolver->data().solution()->valuesGlobal(), 	PETSC_VIEWER_STDOUT_SELF);
   // LOG(DEBUG) << "system matrix of solver: " << *MultiDomainSolver->dataImplicit().systemMatrix();
   // MatView(MultiDomainSolver->dataImplicit().systemMatrix()->valuesGlobal(), PETSC_VIEWER_STDOUT_SELF);


   // run solver
   MultiDomainSolver->advanceTimeSpan();

   // Debug Options
   // VecView(MultiDomainSolver->data().solution()->valuesGlobal(), 	PETSC_VIEWER_STDOUT_SELF);
   // PetscRealView(6, u->values, 0);
   // PetscInt blub;
   // VecGetLocalSize(solution->valuesGlobal(), &blub);
   // LOG(DEBUG) << blub;

   // Get values for braid which were calculated by MultiDomainSolver
   // ierr = VecGetValues(solution->valuesGlobal(), u->size, solution->functionSpace()->meshPartition()->dofNosLocal().data(), u->values); CHKERRQ(ierr);
//    PetscScalar ** val;
//    LOG(DEBUG)<<"test";
//    ierr = VecGetArray(solution->valuesGlobal(), val); CHKERRQ(ierr);
//    LOG(DEBUG)<<"test1";
//    u->values = *val;
//    LOG(DEBUG)<<"test2";
//    ierr = VecRestoreArray(solution->valuesGlobal(), val); CHKERRQ(ierr);
// LOG(DEBUG)<<"test3";
   MultiDomainSolver->getSolution(u->values);

   deltaT = tstop - tstart;
   deltaX = (app->xstop - app->xstart) / (ustop->size - 1.0);

   LOG(DEBUG) << "solved for t: [" << tstart << "," << tstop << "], deltaT: " << deltaT << ", deltaX: " << deltaX << ", size: " << u->size << " solution: " << &(u->values);

   /* Store info on space-time grids visited during the simulation */
   (app->sc_info)[ (2*level) ] = deltaX;
   (app->sc_info)[ (2*level) + 1] = deltaT;

   /* no refinement */
   braid_StepStatusSetRFactor(status, 1);

   return 0;
}

int
my_Init_MD(braid_App     app,
        double        t,
        braid_Vector *u_ptr)
{
   my_Vector *u;
   int nspace = (app->nspace);
   //int    i;
   // double deltaX = (app->xstop - app->xstart) / (nspace - 1.0);
   create_vector_MD(&u, nspace);
   //for(i=0; i < nspace; i++)
   //    {
   //       (u->values)[i] = ((double)braid_Rand())/braid_RAND_MAX;
   //    }

   //LOG(DEBUG) << "app->nspace" << app->nspace;
   /* Allocate vector */
   

   // Different way to initialize
   // /* Initialize vector (with correct boundary conditions) */
   // if(t == 0.0)
   // {
   //    /* Get the solution at time t=0 */
   //    get_solution(u->values, u->size, 0.0, app->xstart, deltaX);
   // }
   // else
   // {
   //    /* Use random values for u(t>0), this measures asymptotic convergence rate */
   //    for(i=0; i < nspace; i++)
   //    {
   //       (u->values)[i] = ((double)braid_Rand())/braid_RAND_MAX;
   //    }
   // }
   // LOG(DEBUG) << "u->size" << u->size;
   //    LOG(DEBUG) << "Braid Init!";
   // LOG(DEBUG) << "---------------";

   PetscInt solver;
   solver=log2(nspace - 1);
   solver=0;
   assert(solver < (*app->MultiDomainSolvers).size());

   std::shared_ptr<typename _braid_App_struct::NestedSolverMD> MultiDomainSolver = (*app->MultiDomainSolvers)[solver];
   //std::shared_ptr<typename Data::PinTMD<typename _braid_App_struct::NestedSolverMD::FunctionSpace>::ScalarFieldVariableType> solution = MultiDomainSolver->data()->transmembranePotential();
   // //auto solution = MultiDomainSolver->data();
   // PetscErrorCode ierr;
   // ierr = VecGetValues(solution->valuesGlobal(), u->size, solution->functionSpace()->meshPartition()->dofNosLocal().data(), u->values); CHKERRQ(ierr);

   // LOG(DEBUG) << "--------------------------------------------------------------";
   // LOG(DEBUG) << "set initial values: " << *solution;

   //PetscScalar a ; 
   //int istart, iend, iterator;
   //VecGetOwnershipRange(MultiDomainSolver->getSolution(u->values),&istart,&iend);
   //LOG(DEBUG) <<istart << ".." << iend;
   //for (iterator=istart; iterator<iend; iterator++) {VecGetValues(MultiDomainSolver->getSolution(), 1, &iterator, &a); u->values [iterator]=a;}
   //int iterator;
   //for (iterator=0; iterator<u->size; iterator++) {std::cout << u->values[iterator];}
   MultiDomainSolver->getSolution(u->values);
   //for (iterator=0; iterator<u->size; iterator++) {std::cout << "lul" << u->values[iterator];}
   *u_ptr = u;

   return 0;
}
