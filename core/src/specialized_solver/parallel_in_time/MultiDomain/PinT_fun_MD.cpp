// most important parallel-in-time structures and functions for XBraid
// more functions in PinT_lib.h

#include "specialized_solver/parallel_in_time/MultiDomain/PinT_fun_MD.h"

#include "specialized_solver/parallel_in_time/MultiDomain/PinT_lib_MD.h"

#include "data_management/specialized_solver/PinT_MD.h"

#include "specialized_solver/parallel_in_time/MultiDomain/PinT_MD_Braid.h"


/* create and allocate a vector */
void
create_vector_MD(my_Vector **u,
              int size)
{
   (*u) = (my_Vector *) malloc(sizeof(my_Vector));
   ((*u)->size)   = size;
   ((*u)->values) = (double *) malloc(size*sizeof(double));
}

int
my_Clone_MD(braid_App     app,
         braid_Vector  u,
         braid_Vector *v_ptr)
{
   my_Vector *v;
   int size = (u->size);
   int i;

   create_vector_MD(&v, size);
   for (i = 0; i < size; i++)
   {
      (v->values)[i] = (u->values)[i];
   }
   *v_ptr = v;

   return 0;
}

int
my_Free_MD(braid_App    app,
        braid_Vector u)
{
   free(u->values);
   free(u);

   return 0;
}

int
my_Sum_MD(braid_App     app,
       double        alpha,
       braid_Vector  x,
       double        beta,
       braid_Vector  y)
{
   int i;
   int size = (y->size);

   for (i = 0; i < size; i++)
   {
      (y->values)[i] = alpha*(x->values)[i] + beta*(y->values)[i];
   }

   return 0;
}


int
my_SpatialNorm_MD(braid_App     app,
               braid_Vector  u,
               double       *norm_ptr)
{
   int    i;
   int size   = (u->size);
   double dot = 0.0;

   for (i = 0; i < size; i++)
   {
      dot += (u->values)[i]*(u->values)[i];
   }
   *norm_ptr = sqrt(dot);

   return 0;
}

int
my_Access_MD(braid_App          app,
          braid_Vector       u,
          braid_AccessStatus astatus)
{
   // int        index, rank, level, done;
   int rank;
   //char       filename[255];
   int        index, level, done;
   double     t, error;

   braid_AccessStatusGetT(astatus, &t);
   braid_AccessStatusGetTIndex(astatus, &index);
   braid_AccessStatusGetLevel(astatus, &level);
   braid_AccessStatusGetDone(astatus, &done);

   /* Print solution to file if simulation is over */
   if(done)
   {
      MPI_Comm_rank( (app->comm), &rank);
      //sprintf(filename, "%s.%07d.%05d", "PinT_diffusion.out", index, rank);
      //save_solution(filename, u->values, u->size, app->xstart,
      //      app->xstop, app->ntime, app->tstart, app->tstop);
      std::shared_ptr<typename _braid_App_struct::NestedSolverMD> MultiDomainSolver = (*app->MultiDomainSolvers)[0];

      //if (rank == 0){
      MultiDomainSolver->setSolution(u->values);
      MultiDomainSolver->printSolution(index, t);
      //}
          
      //for(int i = 0; i < u->size; i++)
      //{
      //   std::cout << "\n" << u->values[i];
      //} 
   }

   /* IF on the finest level AND print_level is high enough AND at the final time,
    * THEN print out the discretization error */
   if( (level == 0) && ((app->print_level) > 0) && (index == app->ntime) )
   {
      error = compute_error_norm(u->values, app->xstart, app->xstop, u->size, t);
      printf("  Discretization error at final time:  %1.4e\n", error);
      fflush(stdout);
   }

   return 0;
}

int
my_BufSize_MD(braid_App           app,
           int                 *size_ptr,
           braid_BufferStatus  bstatus)
{
   int size = (app->nspace);
   *size_ptr = (size+1)*sizeof(double);
   return 0;
}

int
my_BufPack_MD(braid_App           app,
           braid_Vector        u,
           void               *buffer,
           braid_BufferStatus  bstatus)
{
   double *dbuffer = (double *) buffer;
   int i, size = (u->size);

   dbuffer[0] = size;
   for (i = 0; i < size; i++)
   {
      dbuffer[i+1] = (u->values)[i];
   }

   braid_BufferStatusSetSize( bstatus,  (size+1)*sizeof(double));

   return 0;
}

int
my_BufUnpack_MD(braid_App           app,
             void               *buffer,
             braid_Vector       *u_ptr,
             braid_BufferStatus  bstatus)
{
   my_Vector *u = NULL;
   double    *dbuffer = (double *) buffer;
   int        i, size;

   size = dbuffer[0];
   create_vector_MD(&u, size);

   for (i = 0; i < size; i++)
   {
      (u->values)[i] = dbuffer[i+1];
   }
   *u_ptr = u;

   return 0;
}

int
my_Residual_MD(braid_App        app,
            braid_Vector     ustop,
            braid_Vector     r,
            braid_StepStatus status)
{
   // double tstart;             /* current time */
   // double tstop;              /* evolve to this time*/
   // int i;
   // double x, deltaX, deltaT;

   // braid_StepStatusGetTstartTstop(status, &tstart, &tstop);
   // deltaT = tstop - tstart;
   // deltaX = (app->xstop - app->xstart) / (ustop->size - 1.0);

   // /* Set up matrix stencil for 1D heat equation*/
   // compute_stencil(deltaX, deltaT, app->matrix);

   // /* Residual r = A*xstop - r - forcing - boundary
   //  *   note: there are no boundary terms here */
   // matvec_tridiag(ustop->values, app->g, ustop->size, app->matrix);
   // x = app->xstart;
   // for(i = 0; i < r->size; i++)
   // {
   //    r->values[i] = app->g[i] - r->values[i] - deltaT*forcing(tstop, x);
   //    x = x + deltaX;
   // }

   return 0;
}

/* Bilinear Coarsening */
int
my_Coarsen_MD(braid_App              app,
           braid_Vector           fu,
           braid_Vector          *cu_ptr,
           braid_CoarsenRefStatus status)
{

   int csize, level;
   my_Vector *v;

   /* This returns the level for fu */
   braid_CoarsenRefStatusGetLevel(status, &level);

   /* Stop spatial coarsening/interpolation after reaching a grid of 3 points.
    * This is the smallest spatial grid (2 boundary points, one true DOF). */
   if( level < floor(log2(app->nspace)) - 1 )
   {
      csize = ((fu->size) - 1)/2 + 1;
      create_vector_MD(&v, csize);
      coarsen_1D(v->values, fu->values, csize, fu->size);
   }
   else
   {
      /* No coarsening, clone the vector */
      my_Clone_MD(app, fu, &v);
   }

   *cu_ptr = v;

   return 0;
}

/* Bilinear interpolation */
int
my_Interp_MD(braid_App              app,
          braid_Vector           cu,
          braid_Vector          *fu_ptr,
          braid_CoarsenRefStatus status)
{

   int fsize, level;
   my_Vector *v;

   /* This returns the level for fu_ptr */
   braid_CoarsenRefStatusGetLevel(status, &level);

   /* Stop spatial coarsening/interpolation after reaching a grid of 3 points.
    * This is the smallest spatial grid (2 boundary points, one true DOF). */
   if( level < floor(log2(app->nspace)) - 1 )
   {
      fsize = (cu->size - 1)*2 + 1;
      create_vector_MD(&v, fsize);
      interpolate_1D(cu->values, v->values, cu->size, fsize);
   }
   else
   {
      /* No refinement, clone the vector */
      my_Clone_MD(app, cu, &v);
   }

*fu_ptr = v;

   return 0;
}

