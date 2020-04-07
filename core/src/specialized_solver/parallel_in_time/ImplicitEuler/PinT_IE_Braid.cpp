// specific functions for XBraid with Implicit Euler
// defines my_Step (integration step)

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
   std::shared_ptr<typename _braid_App_struct::NestedSolver> implicitEulerSolver = (*app->implicitEulerSolvers)[solver];
   std::shared_ptr<typename Data::PinTIE<typename _braid_App_struct::NestedSolver::FunctionSpace>::ScalarFieldVariableType> solution = implicitEulerSolver->data().solution();
   assert(u->size == solution->nDofsGlobal());

   // Set the initial guess for the solver
   PetscErrorCode ierr;
   ierr = VecSetValues(solution->valuesGlobal(), u->size, solution->functionSpace()->meshPartition()->dofNosLocal().data(), u->values, INSERT_VALUES); CHKERRQ(ierr);
   VecAssemblyBegin(solution->valuesGlobal());
   VecAssemblyEnd(solution->valuesGlobal());

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

   // Get values for braid which were calculated by implicitEulerSolver
   ierr = VecGetValues(solution->valuesGlobal(), u->size, solution->functionSpace()->meshPartition()->dofNosLocal().data(), u->values); CHKERRQ(ierr);
//    PetscScalar ** val;
//    LOG(DEBUG)<<"test";
//    ierr = VecGetArray(solution->valuesGlobal(), val); CHKERRQ(ierr);
//    LOG(DEBUG)<<"test1";
//    u->values = *val;
//    LOG(DEBUG)<<"test2";
//    ierr = VecRestoreArray(solution->valuesGlobal(), val); CHKERRQ(ierr);
// LOG(DEBUG)<<"test3";
   deltaT = tstop - tstart;
   deltaX = (app->xstop - app->xstart) / (ustop->size - 1.0);

   LOG(DEBUG) << "solved for t: [" << tstart << "," << tstop << "], deltaT: " << deltaT << ", deltaX: " << deltaX << ", size: " << u->size << " solution: " << *solution;

   /* Store info on space-time grids visited during the simulation */
   (app->sc_info)[ (2*level) ] = deltaX;
   (app->sc_info)[ (2*level) + 1] = deltaT;

   /* no refinement */
   braid_StepStatusSetRFactor(status, 1);

   return 0;
}
