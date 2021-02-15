// most important parallel-in-time structures and functions for XBraid
// more functions in PinT_lib.h

#include "specialized_solver/parallel_in_time/ImplicitEuler/PinT_fun_IE.h"

#include "specialized_solver/parallel_in_time/ImplicitEuler/PinT_lib_IE.h"

// /* create and allocate a vector */
// void
// create_vector(my_Vector **u,
//               int size)
// {
//    (*u) = (my_Vector *) malloc(sizeof(my_Vector));
//    ((*u)->size)   = size;
//    ((*u)->values) = (double *) malloc(size*sizeof(double));
// }

// int
// my_Clone(braid_App     app,
//          braid_Vector  u,
//          braid_Vector *v_ptr)
// {
//    my_Vector *v;
//    int size = (u->size);
//    int i;

//    create_vector(&v, size);
//    for (i = 0; i < size; i++)
//    {
//       (v->values)[i] = (u->values)[i];
//    }
//    *v_ptr = v;

//    return 0;
// }

// int
// my_Free(braid_App    app,
//         braid_Vector u)
// {
//    free(u->values);
//    free(u);

//    return 0;
// }

// int
// my_Sum(braid_App     app,
//        double        alpha,
//        braid_Vector  x,
//        double        beta,
//        braid_Vector  y)
// {
//    int i;
//    int size = (y->size);

//    for (i = 0; i < size; i++)
//    {
//       (y->values)[i] = alpha*(x->values)[i] + beta*(y->values)[i];
//    }

//    return 0;
// }


// int
// my_SpatialNorm(braid_App     app,
//                braid_Vector  u,
//                double       *norm_ptr)
// {
//    int    i;
//    int size   = (u->size);
//    double dot = 0.0;

//    for (i = 0; i < size; i++)
//    {
//       dot += (u->values)[i]*(u->values)[i];
//    }
//    *norm_ptr = sqrt(dot);

//    return 0;
// }

// int
// my_Access(braid_App          app,
//           braid_Vector       u,
//           braid_AccessStatus astatus)
// {
//    // int        index, rank, level, done;
//    int rank;
//    char       filename[255];
//    int        index, level, done;
//    double     t, error;

//    braid_AccessStatusGetT(astatus, &t);
//    braid_AccessStatusGetTIndex(astatus, &index);
//    braid_AccessStatusGetLevel(astatus, &level);
//    braid_AccessStatusGetDone(astatus, &done);

//    /* Print solution to file if simulation is over */
//    if(done)
//    {
//       //MPI_Comm_rank( (app->comm), &rank);
//       //sprintf(filename, "%s.%07d.%05d", "PinT_diffusion.out", index, rank);
//       //save_solution(filename, u->values, u->size, app->xstart,
//       //      app->xstop, app->ntime, app->tstart, app->tstop);
//    }

//    /* IF on the finest level AND print_level is high enough AND at the final time,
//     * THEN print out the discretization error */
//    if( (level == 0) && ((app->print_level) > 0) && (index == app->ntime) )
//    {
//       error = compute_error_norm(u->values, app->xstart, app->xstop, u->size, t);
//       printf("  Discretization error at final time:  %1.4e\n", error);
//       fflush(stdout);
//    }

//    return 0;
// }

// int
// my_BufSize(braid_App           app,
//            int                 *size_ptr,
//            braid_BufferStatus  bstatus)
// {
//    int size = (app->nspace);
//    *size_ptr = (size+1)*sizeof(double);
//    return 0;
// }

// int
// my_BufPack(braid_App           app,
//            braid_Vector        u,
//            void               *buffer,
//            braid_BufferStatus  bstatus)
// {
//    double *dbuffer = (double *) buffer;
//    int i, size = (u->size);

//    dbuffer[0] = size;
//    for (i = 0; i < size; i++)
//    {
//       dbuffer[i+1] = (u->values)[i];
//    }

//    braid_BufferStatusSetSize( bstatus,  (size+1)*sizeof(double));

//    return 0;
// }

// int
// my_BufUnpack(braid_App           app,
//              void               *buffer,
//              braid_Vector       *u_ptr,
//              braid_BufferStatus  bstatus)
// {
//    my_Vector *u = NULL;
//    double    *dbuffer = (double *) buffer;
//    int        i, size;

//    size = dbuffer[0];
//    create_vector(&u, size);

//    for (i = 0; i < size; i++)
//    {
//       (u->values)[i] = dbuffer[i+1];
//    }
//    *u_ptr = u;

//    return 0;
// }

// int
// my_Residual(braid_App        app,
//             braid_Vector     ustop,
//             braid_Vector     r,
//             braid_StepStatus status)
// {
//    double tstart;             /* current time */
//    double tstop;              /* evolve to this time*/
//    int i;
//    double x, deltaX, deltaT;

//    braid_StepStatusGetTstartTstop(status, &tstart, &tstop);
//    deltaT = tstop - tstart;
//    deltaX = (app->xstop - app->xstart) / (ustop->size - 1.0);

//    /* Set up matrix stencil for 1D heat equation*/
//    compute_stencil(deltaX, deltaT, app->matrix);

//    /* Residual r = A*xstop - r - forcing - boundary
//     *   note: there are no boundary terms here */
//    matvec_tridiag(ustop->values, app->g, ustop->size, app->matrix);
//    x = app->xstart;
//    for(i = 0; i < r->size; i++)
//    {
//       r->values[i] = app->g[i] - r->values[i] - deltaT*forcing(tstop, x);
//       x = x + deltaX;
//    }

//    return 0;
// }

// /* Bilinear Coarsening */
// int
// my_Coarsen(braid_App              app,
//            braid_Vector           fu,
//            braid_Vector          *cu_ptr,
//            braid_CoarsenRefStatus status)
// {

//    int csize, level;
//    my_Vector *v;

//    /* This returns the level for fu */
//    braid_CoarsenRefStatusGetLevel(status, &level);

//    /* Stop spatial coarsening/interpolation after reaching a grid of 3 points.
//     * This is the smallest spatial grid (2 boundary points, one true DOF). */
//    if( level < floor(log2(app->nspace)) - 1 )
//    {
//       csize = ((fu->size) - 1)/2 + 1;
//       create_vector(&v, csize);
//       coarsen_1D(v->values, fu->values, csize, fu->size);
//    }
//    else
//    {
//       /* No coarsening, clone the vector */
//       my_Clone(app, fu, &v);
//    }

//    *cu_ptr = v;

//    return 0;
// }

// /* Bilinear interpolation */
// int
// my_Interp(braid_App              app,
//           braid_Vector           cu,
//           braid_Vector          *fu_ptr,
//           braid_CoarsenRefStatus status)
// {

//    int fsize, level;
//    my_Vector *v;

//    /* This returns the level for fu_ptr */
//    braid_CoarsenRefStatusGetLevel(status, &level);

//    /* Stop spatial coarsening/interpolation after reaching a grid of 3 points.
//     * This is the smallest spatial grid (2 boundary points, one true DOF). */
//    if( level < floor(log2(app->nspace)) - 1 )
//    {
//       fsize = (cu->size - 1)*2 + 1;
//       create_vector(&v, fsize);
//       interpolate_1D(cu->values, v->values, cu->size, fsize);
//    }
//    else
//    {
//       /* No refinement, clone the vector */
//       my_Clone(app, cu, &v);
//    }

//    *fu_ptr = v;

//    return 0;
// }

/* create and allocate a vector */
void
create_vector(my_Vector **u,
              int size, braid_App app)
{
   (*u) = (my_Vector *) malloc(sizeof(my_Vector));
 
   VecCreate(app->comm,&(*u)->values);
   VecSetSizes((*u)->values, size, PETSC_DECIDE);
   VecSetFromOptions((*u)->values);
 //((*u)->size)   = size;
 //((*u)->values) = (double *) malloc(size*sizeof(double));
}

int
my_Clone(braid_App     app,
         braid_Vector  u,
         braid_Vector *v_ptr)
{
   my_Vector *v;
   //DM dmf;
   //double *_b;
   (v) = (my_Vector *) malloc(sizeof(my_Vector));
   //int i;
   //int size, sizeu;
   //VecGetOwnershipRange(u->values, &size, &sizeu);
   //VecGetSize(u->values, &size);
   //LOG(DEBUG) << size;
   //VecGetLocalSize(u->values, &sizeu);
   //LOG(DEBUG) << sizeu;
   //VecView(u->values, PETSC_VIEWER_STDOUT_WORLD);
   //VecGetDM(u->values, &dmf);
   //VecGetArray(u->values, &_b);
   //for (i = 0; i < size; i++)
   // {
   //    LOG(DEBUG) << i << " " << _b[i];
   // }
   
   //DMCreateGlobalVector(dmf, &v->values);
   VecDuplicate(u->values, &v->values);
   VecCopy(u->values,v->values);
   *v_ptr = v;

   return 0;
}

int
my_Free(braid_App    app,
        braid_Vector u)
{
   VecDestroy(&u->values);
   free(u);

   return 0;
}

int
my_Sum(braid_App     app,
       double        alpha,
       braid_Vector  x,
       double        beta,
       braid_Vector  y)
{
   Vec  help;
   VecDuplicate(x->values, &help);
   VecCopy(x->values, help);
   VecAXPBY(y->values, alpha, beta, help);
   return 0;
}


int
my_SpatialNorm(braid_App     app,
               braid_Vector  u,
               double       *norm_ptr)
{
   VecNorm(u->values, NORM_2, norm_ptr);
   return 0;
}

int
my_Access(braid_App          app,
          braid_Vector       u,
          braid_AccessStatus astatus)
{
   // int        index, rank, level, done;
   int rank,size;
   char       filename[255];
   int        index, level, done;
   double     t, error;

   braid_AccessStatusGetT(astatus, &t);
   braid_AccessStatusGetTIndex(astatus, &index);
   braid_AccessStatusGetLevel(astatus, &level);
   braid_AccessStatusGetDone(astatus, &done);

   /* Print solution to file if simulation is over */
   if(done)
   {
      //MPI_Comm_rank( (app->comm), &rank);
      //sprintf(filename, "%s.%07d.%05d", "PinT_diffusion.out", index, rank);
      //PetscScalar    *array;
      //VecGetSize(u->values, &size);
      //VecGetArray(u->values, &array);
      //save_solution(filename, array, size, app->xstart,
      //      app->xstop, app->ntime, app->tstart, app->tstop);
      //VecRestoreArray(u->values, &array);
   }

   /* IF on the finest level AND print_level is high enough AND at the final time,
    * THEN print out the discretization error */
   if( (level == 0) && ((app->print_level) > 0) && (index == app->ntime) )
   {
      //error = compute_error_norm(u->values, app->xstart, app->xstop, u->size, t);
      //printf("  Discretization error at final time:  %1.4e\n", error);
      //fflush(stdout);
   }

   return 0;
}

int
my_BufSize(braid_App           app,
           int                 *size_ptr,
           braid_BufferStatus  bstatus)
{
   int size = (app->nspace);
   *size_ptr = (size+1)*sizeof(double);
   return 0;
}

int
my_BufPack(braid_App           app,
           braid_Vector        u,
           void               *buffer,
           braid_BufferStatus  bstatus)
{
   double *dbuffer = (double *) buffer;
   int i, size, low, high;
   VecGetOwnershipRange(u->values, &low, &high);
   VecGetLocalSize(u->values, &size);
   dbuffer[0] = size;
   double    *_a;
   VecGetArray(u->values,&_a);
   //for (i = low; i < high; i++) dbuffer[i+1-low] = _a[i];
   for (i = low; i < high; i++) dbuffer[i+1-low] = _a[i-low];
   VecRestoreArray(u->values,&_a);

   braid_BufferStatusSetSize( bstatus,  (size+1)*sizeof(double));

   return 0;
}

int
my_BufUnpack(braid_App           app,
             void               *buffer,
             braid_Vector       *u_ptr,
             braid_BufferStatus  bstatus)
{
   my_Vector *u = NULL;
   //(u) = (my_Vector *) malloc(sizeof(my_Vector));
   //VecDuplicate(u_ptr->values, &u->values);
   double    *dbuffer = (double *) buffer;
   int        iterator, size, istarti, iendi;

   size = dbuffer[0];
   create_vector(&u, size, app);

   VecGetOwnershipRange(u->values,&istarti,&iendi);
   //for (iterator=0; iterator<u->size; iterator++) {std::cout<<u->values[iterator] << "/n";}
   for (iterator=istarti; iterator<iendi; iterator++) {VecSetValues(u->values, 1, &iterator, &dbuffer[iterator+1], INSERT_VALUES);}
   //MPI_Barrier(app->comm);
   VecAssemblyBegin(u->values);
   VecAssemblyEnd(u->values);

   *u_ptr = u;

   return 0;
}

// int
// my_Residual(braid_App        app,
//             braid_Vector     ustop,
//             braid_Vector     r,
//             braid_StepStatus status)
// {
//    double tstart;             /* current time */
//    double tstop;              /* evolve to this time*/
//    int i;
//    double x, deltaX, deltaT;

//    braid_StepStatusGetTstartTstop(status, &tstart, &tstop);
//    deltaT = tstop - tstart;
//    deltaX = (app->xstop - app->xstart) / (ustop->size - 1.0);

//    /* Set up matrix stencil for 1D heat equation*/
//    compute_stencil(deltaX, deltaT, app->matrix);

//    /* Residual r = A*xstop - r - forcing - boundary
//     *   note: there are no boundary terms here */
//    matvec_tridiag(ustop->values, app->g, ustop->size, app->matrix);
//    x = app->xstart;
//    for(i = 0; i < r->size; i++)
//    {
//       r->values[i] = app->g[i] - r->values[i] - deltaT*forcing(tstop, x);
//       x = x + deltaX;
//    }

//    return 0;
// }

/* Bilinear Coarsening */
int
my_Coarsen(braid_App              app,
           braid_Vector           fu,
           braid_Vector          *cu_ptr,
           braid_CoarsenRefStatus status)
{

   int level;
   PetscErrorCode ierr;
   my_Vector *v;
   v = (my_Vector *) malloc(sizeof(my_Vector));
   //Vec v->values;
   //VecCreate(app->comm, &v->values);
   //VecSetType(v->values, VECMPI);
   //VecAssemblyBegin(H);
   //VecAssemblyEnd(H);

   /* This returns the level for fu */
   braid_CoarsenRefStatusGetLevel(status, &level);

   /* Stop spatial coarsening/interpolation after reaching a grid of 3 points.
    * This is the smallest spatial grid (2 boundary points, one true DOF). */
   if( level < floor(log2(app->nspace)) - 1 )
   {
      int csize;
      //v = (my_Vector *) malloc(sizeof(my_Vector));
      //DM dmf, fdm;
      int size_fu, size_v, i;
      Mat A;
      //double    *_fu;
      //double    *_v;
      //VecGetDM(fu->values, &fdm);
      //DMCoarsen(fdm, app->comm, &dmf);
      //DMCreateGlobalVector(dmf, &v->values);
      VecGetSize(fu->values, &size_fu);
      //VecGetSize(v->values, &size_v);
      //VecGetArray(fu->values, &_fu);
      // VECGET ARRAY GIVES LOCAL NOT GLOBAL!!!
      //VecGetArray(v->values, &_v);
      //for (i = 0; i < size_fu; i++)
      //{
      // LOG(DEBUG) << i << " asdsad " << _fu[i];
      //}
      //coarsen_1D(_v, _fu, size_v, size_fu);
      //VecRestoreArray(fu->values, &_fu);
      //VecRestoreArray(v->values, &_v);
      //VecView(v->values, PETSC_VIEWER_STDOUT_WORLD);
      //VecView(fu->values, PETSC_VIEWER_STDOUT_WORLD);
      PetscScalar k;
      k=1;
      PetscInt    vi[3];
      PetscScalar vs[3];
      vs[0]=0.25;
      vs[1]=0.5;
      vs[2]=0.25;
      csize = (size_fu - 1)/2 + 1;
      //VecSetSizes(v->values, PETSC_DECIDE, csize);
      MatCreateAIJ(app->comm,PETSC_DECIDE,PETSC_DECIDE,csize,size_fu,3,NULL,1,NULL,&A);
      MatSetUp(A);
      for (i = 0; i < csize; i++)
      {
         if (i==0)
         {
            MatSetValues(A, 1, &i, 1, &i, &k, INSERT_VALUES);
         }
         else if (i==(csize-1))
         {
            int a;
            a=size_fu-1;
            MatSetValues(A, 1, &i, 1, &a, &k, INSERT_VALUES);
         }
         else
         {
            int j;
            j=i*2;
            vi[0]=j-1;
            vi[1]=j;
            vi[2]=j+1;
            MatSetValues(A, 1, &i, 3, vi, vs, INSERT_VALUES);
         }
      }
      MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
      MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
      MatCreateVecs(A, NULL, &v->values);
      //MatView(A, PETSC_VIEWER_STDOUT_WORLD);
      //VecView(fu->values, PETSC_VIEWER_STDOUT_WORLD);
      //PetscInt row,col;
      //MatGetOwnershipRange(A,&row,&col);
      //LOG(DEBUG)<< row << "   " << col;
      ierr = MatMult(A, fu->values, v->values); CHKERRQ(ierr);
      //VecView(v->values, PETSC_VIEWER_STDOUT_WORLD);
   }
   else
   {
      /* No coarsening, clone the vector */
      //v = (my_Vector *) malloc(sizeof(my_Vector));
      my_Clone(app, fu, &v);
   }

   *cu_ptr = v;

   return 0;
}

/* Bilinear interpolation */
int
my_Interp(braid_App              app,
          braid_Vector           cu,
          braid_Vector          *fu_ptr,
          braid_CoarsenRefStatus status)
{

   int level;
   int size_cu, size_v, i;
   my_Vector *v;
   v = (my_Vector *) malloc(sizeof(my_Vector));
   VecCreate(app->comm, &v->values);
   VecSetType(v->values, VECMPI);
   //DM fdm, cdm;
   //double    *_cu;
   //double    *_v;

   /* This returns the level for fu_ptr */
   braid_CoarsenRefStatusGetLevel(status, &level);

   /* Stop spatial coarsening/interpolation after reaching a grid of 3 points.
    * This is the smallest spatial grid (2 boundary points, one true DOF). */
   if( level < floor(log2(app->nspace)) - 1 )
   {
      //v = (my_Vector *) malloc(sizeof(my_Vector));
      //VecGetDM(cu->values, &cdm);
      //DMRefine(cdm, app->comm, &fdm);
      //DMCreateGlobalVector(fdm, &v->values);
      //VecGetSize(cu->values, &size_cu);
      //VecGetSize(v->values, &size_v);
      //VecGetArray(cu->values, &_cu);
      //VecGetArray(v->values, &_v);
      //interpolate_1D(_cu, _v, size_cu, size_v);
      //VecRestoreArray(cu->values, &_cu);
      //VecRestoreArray(v->values, &_v);

      VecGetSize(cu->values, &size_cu);
      Mat A;
      PetscScalar k;
      k=1;
      PetscInt    vi[2];
      PetscScalar vs[2];
      vs[0]=0.5;
      vs[1]=0.5;
      size_v = (size_cu - 1)*2 + 1;
      VecSetSizes(v->values, PETSC_DECIDE, size_v);
      MatCreateAIJ(app->comm,PETSC_DECIDE,PETSC_DECIDE,size_v,size_cu,3,NULL,1,NULL,&A);
      MatSetUp(A);
      for (i = 0; i < size_v; i++)
      {
         if (i==0)
         {
            MatSetValues(A, 1, &i, 1, &i, &k, INSERT_VALUES);
         }
         else if (i==(size_v-1))
         {
            int a;
            a=size_cu-1;
            MatSetValues(A, 1, &i, 1, &a, &k, INSERT_VALUES);
         }
         else
         {
            int j,l;
            j=i/2;
            l=(i+1)/2;
            if(i%2 == 1)
            {
               MatSetValues(A, 1, &i, 1, &j, &k, INSERT_VALUES);
            }
            else
            {
               vi[0]=j;
               vi[1]=l;
               MatSetValues(A, 1, &i, 2, vi, vs, INSERT_VALUES);
            }
         }
      }
      MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
      MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
      //MatView(A, PETSC_VIEWER_STDOUT_WORLD);
      //VecView(cu->values, PETSC_VIEWER_STDOUT_WORLD);
      MatMult(A, cu->values, v->values);
      //VecView(v->values, PETSC_VIEWER_STDOUT_WORLD);
   }
   else
   {
      /* No refinement, clone the vector */
      //v = (my_Vector *) malloc(sizeof(my_Vector));
      my_Clone(app, cu, &v);
   }

   *fu_ptr = v;

   return 0;
}