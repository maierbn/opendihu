// specific functions for XBraid with Implicit Euler
// defines my_Step (integration step)
#include <mpi.h>

#include "specialized_solver/parallel_in_time/ImplicitEuler/PinT_IE_Braid.h"

#include "specialized_solver/parallel_in_time/PinT_fun.h"
#include "specialized_solver/parallel_in_time/PinT_lib.h"
#include "data_management/specialized_solver/PinT_IE.h"

int my_Step(braid_App        app,
            braid_Vector     ustop,
            braid_Vector     fstop,
            braid_Vector     u,
            braid_StepStatus status)
{
   PetscReal tstart;             /* current time */
   PetscReal tstop;              /* evolve to this time*/
   PetscInt level, i, solver;
   PetscReal deltaX, deltaT;
   // PetscReal * help;

   // get level and start and stop time from braid status
   braid_StepStatusGetLevel(status, &level);
   braid_StepStatusGetTstartTstop(status, &tstart, &tstop);

   // determine, which solver is needed (depends on size)
   solver=log2(u->size - 1);
   assert(solver < (*app->implicitEulerSolvers).size());

   /* XBraid forcing */

   if(fstop != NULL)
   {
      for(i = 0; i < u->size; i++)
      {
         u->values[i] = u->values[i] + fstop->values[i];
      }
   }

   // get alias variables
   std::shared_ptr<typename _braid_App_struct::NestedSolverIE> implicitEulerSolver = (*app->implicitEulerSolvers)[solver];
   std::shared_ptr<typename Data::PinTIE<typename _braid_App_struct::NestedSolverIE::FunctionSpace>::ScalarFieldVariableType> solution = implicitEulerSolver->data().solution();

   assert(u->size == solution->nDofsGlobal());

   // Set the initial guess for the solver
   PetscErrorCode ierr;
   //ierr = VecSetValues(solution->valuesGlobal(), u->size, solution->functionSpace()->meshPartition()->dofNosLocal().data(), u->values, INSERT_VALUES); CHKERRQ(ierr);
   //VecAssemblyBegin(solution->valuesGlobal());
   //VecAssemblyEnd(solution->valuesGlobal());

   //VecView(implicitEulerSolver->data().solution()->valuesGlobal(), PETSC_VIEWER_STDOUT_WORLD);

   PetscInt iterator;
   int istarti, iendi, itotal;
   itotal = iendi-istarti;
   VecGetOwnershipRange(solution->valuesGlobal(),&istarti,&iendi);
   //for (iterator=0; iterator<u->size; iterator++) {std::cout<<u->values[iterator] << "/n";}
   for (iterator=istarti; iterator<iendi; iterator++) {VecSetValues(solution->valuesGlobal(), 1, &iterator, &u->values[iterator], INSERT_VALUES);}
   //MPI_Barrier(app->comm);
   VecAssemblyBegin(solution->valuesGlobal());
   VecAssemblyEnd(solution->valuesGlobal());

   //VecView(implicitEulerSolver->data().solution()->valuesGlobal(), PETSC_VIEWER_STDOUT_WORLD);
   //int rank, size;
   //size=0;
   //MPI_Comm_size(app->comm, &size);
   //if (size > 1)
  // {
   //   MPI_Comm_rank(app->comm, &rank);
   //   MPI_Bcast(u->values, itotal, MPI_DOUBLE, rank, app->comm);
   //   MPI_Barrier(app->comm);
   //}



   // set time span for the solver, which is calculated by braid status
   implicitEulerSolver->setTimeSpan(tstart, tstop);
   implicitEulerSolver->setNumberTimeSteps(1);
   implicitEulerSolver->setSystemMatrix((tstop-tstart)/1);

   // Debug Options
   LOG(DEBUG) << "--------------------------------------------------------------";
   LOG(DEBUG) << "level: " << level << ", solver: " << solver << ", size: " << u->size << ", t: [" << tstart << "," << tstop << "], before implicit euler:" << *solution;

   // PetscRealView(u->size, u->values, 0);
   // VecView(implicitEulerSolver->data().solution()->valuesGlobal(), 	PETSC_VIEWER_STDOUT_SELF);
   // LOG(DEBUG) << "system matrix of solver: " << *implicitEulerSolver->dataImplicit().systemMatrix();
   // MatView(implicitEulerSolver->dataImplicit().systemMatrix()->valuesGlobal(), PETSC_VIEWER_STDOUT_SELF);


   // run solver
   implicitEulerSolver->advanceTimeSpan();

   // Debug Options
   // VecView(implicitEulerSolver->data().solution()->valuesGlobal(), 	PETSC_VIEWER_STDOUT_SELF);
   // PetscRealView(6, u->values, 0);
   // PetscInt blub;
   // VecGetLocalSize(solution->valuesGlobal(), &blub);
   // LOG(DEBUG) << blub;

   //VecScatter ctx;
   //Vec vout;
   //if (app->testscatter == 0){
   //   app->testscatter=u->size;
   //   VecScatterCreateToAll(solution->valuesGlobal(),&app->vecscatter,&vout);
   //
   //}
   //else if (app->testscatter!=u->size){
   //   app->testscatter=u->size;
   //   VecScatterDestroy(&app->vecscatter);
   //   VecScatterCreateToAll(solution->valuesGlobal(),&app->vecscatter,&vout);
   //
   //}
   //VecScatterCreateToAll(solution->valuesGlobal(),&app->vecscatter,&vout);
   // scatter as many times as you need
   //VecScatterBegin(app->vecscatter,solution->valuesGlobal(),vout,INSERT_VALUES,SCATTER_FORWARD);
   //VecScatterEnd(app->vecscatter,solution->valuesGlobal(),vout,INSERT_VALUES,SCATTER_FORWARD);
   // destroy scatter context and local vector when no longer needed
   //double    *_a;
   //VecGetArray(vout,&_a);
   //for (i = 0; i < u->size; i++) u->values[i] = _a[i];
   //VecRestoreArray(vout,&_a);
   //VecScatterDestroy(&app->vecscatter);
   //VecDestroy(&vout);

   //int istart, iend;
   //VecGetOwnershipRange(solution->valuesGlobal(),&istart,&iend);

   // Get values for braid which were calculated by implicitEulerSolver
   //PetscScalar arr;
   //for (iterator=istart; iterator<iend; iterator++) {VecGetValues(solution->valuesGlobal(), 1, &iterator, &arr); u->values[iterator]=arr;}
   //MPI_Barrier(app->comm);
   ierr = VecGetValues(solution->valuesGlobal(), u->size, solution->functionSpace()->meshPartition()->dofNosLocal().data(), u->values); CHKERRQ(ierr);
   //std::cout << "asd \n";
   //for (iterator=0; iterator<u->size; iterator++) {std::cout << "debug"; std::cout<<u->values[iterator] << "/n";}
   //MPI_Allreduce(MPI_IN_PLACE, u->values, u->size, MPI_DOUBLE, MPI_MAX, app->comm);
   //MPI_Allgather(u->values, )
   //for (iterator=istart; iterator<iend; iterator++) {std::cout << "| \n" << u->values[iterator];}

   //VecScatterDestroy(&ctx);
   //VecDestroy(&vout);

   //for (iterator=0; iterator<u->size; iterator++) {std::cout<<u->values[iterator] << "/n";}
   deltaT = tstop - tstart;
   deltaX = (app->xstop - app->xstart) / (ustop->size - 1.0);

   LOG(DEBUG) << "solved for t: [" << tstart << "," << tstop << "], deltaT: " << deltaT << ", deltaX: " << deltaX << ", size: " << u->size << " solution: " << *solution;

   //VecView(implicitEulerSolver->data().solution()->valuesGlobal(), PETSC_VIEWER_STDOUT_WORLD);
   /* Store info on space-time grids visited during the simulation */
   (app->sc_info)[ (2*level) ] = deltaX;
   (app->sc_info)[ (2*level) + 1] = deltaT;

   /* no refinement */
   braid_StepStatusSetRFactor(status, 1);

   return 0;
}

int
my_Init(braid_App     app,
        double        t,
        braid_Vector *u_ptr)
{
   my_Vector *u;
   int nspace = (app->nspace);
   // int    i;
   // double deltaX = (app->xstop - app->xstart) / (nspace - 1.0);

   /* Allocate vector */
   create_vector(&u, nspace);

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

   PetscInt solver;
   solver=log2(nspace - 1);
   assert(solver < (*app->implicitEulerSolvers).size());

   std::shared_ptr<typename _braid_App_struct::NestedSolverIE> implicitEulerSolver = (*app->implicitEulerSolvers)[solver];
   std::shared_ptr<typename Data::PinTIE<typename _braid_App_struct::NestedSolverIE::FunctionSpace>::ScalarFieldVariableType> solution = implicitEulerSolver->data().solution();

   PetscScalar a ; 
   int istart, iend, iterator;
   VecGetOwnershipRange(solution->valuesGlobal(),&istart,&iend);
   LOG(DEBUG) <<istart << ".." << iend;
   for (iterator=0; iterator<u->size; iterator++) {u->values[iterator]=0;}
   for (iterator=istart; iterator<iend; iterator++) {VecGetValues(solution->valuesGlobal(), 1, &iterator, &a); u->values [iterator]=a;}
   //for (iterator=istart; iterator<iend; iterator++) {u->values[iterator]=0;}
   //ierr = VecGetValues(solution->valuesGlobal(), u->size, solution->functionSpace()->meshPartition()->dofNosLocal().data(), u->values); CHKERRQ(ierr);

   LOG(DEBUG) << "--------------------------------------------------------------";
   LOG(DEBUG) << "set initial values: " << *solution;

   *u_ptr = u;

   return 0;
}
