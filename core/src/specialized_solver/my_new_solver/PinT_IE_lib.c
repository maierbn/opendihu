/*BHEADER**********************************************************************
 * Copyright (c) 2013, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory. Written by
 * Jacob Schroder, Rob Falgout, Tzanio Kolev, Ulrike Yang, Veselin
 * Dobrev, et al. LLNL-CODE-660355. All rights reserved.
 *
 * This file is part of XBraid. For support, post issues to the XBraid Github page.
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License (as published by the Free Software
 * Foundation) version 2.1 dated February 1999.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE. See the terms and conditions of the GNU General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc., 59
 * Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 ***********************************************************************EHEADER*/


/**
 *  This file contains library functions for ex-02.c and ex-02-serial.c.
 *  Together, these files solve the 1D heat equation.
 *
 *  For more details on the discretization, see the header comment in ex-02.c.
 **/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include <petscdm.h>
#include <petscdmda.h>
#include <petscts.h>
#include <petscdraw.h>

#ifdef M_PI
   #define PI M_PI
#else
   #define PI 3.14159265358979
#endif

/*--------------------------------------------------------------------------
 * Routines defining solution, forcing term and how to compute the error
 *--------------------------------------------------------------------------*/

/* Exact solution */
double exact(PetscReal t, PetscReal x)
{
    return sin(x)*cos(t);
}

/* Initialize array of values to solution at time t */
void
get_solution(PetscReal * values,
             PetscInt      size,
             PetscReal   t,
             PetscReal   xstart,
             PetscReal   deltaX)
{
   PetscInt i;

   for(i=0; i < size; i++)
   {
      values[i] = exact(t, xstart);
      xstart += deltaX;
   }
}

/* Forcing term F(t, x) for PDE,  u_t = u_xx + F(t,x) */
double forcing(PetscReal t, PetscReal x)
{
   return (-1.0)*sin(x)*sin(t) + sin(x)*cos(t);
}

/* Compute L2-norm of the error at a point in time */
double compute_error_norm(double * values,
                          double   xstart,
                          double   xstop,
                          int      nspace,
                          double   t)
{
   int i;
   double deltaX = (xstop - xstart) / (nspace - 1.0);
   double x = xstart;
   double error = 0.0;

   for(i = 0; i < nspace; i++)
   {
      error = error + pow( values[i] - exact(t, x), 2);
      x += deltaX;
   }

   return sqrt( error*deltaX );
}


/*--------------------------------------------------------------------------
 * Routines for computing time-steps and residuals
 *--------------------------------------------------------------------------*/

/* Helper function for tridiagonal solver */
PetscReal dabs(PetscReal x)
{
   if (x < 0.0)
   {
      return -1.0;
   }
   else if (x > 0.0)
   {
      return 1.0;
   }
   else
   {
      return 0.0;
   }
}


/* Helper function for Step: Tridiagonal system solver (Thomas algorithm) */
void
solve_tridiag(PetscReal *x, PetscReal *g, PetscInt N, PetscReal* matrix)
{
   /**
    * solves Ax = v where A is a tridiagonal matrix with stencil [ a, b, c].
    *
    * There is a built in assumption that the first and last rows are the
    * identity (boundary conditions)
    *
    * Input
    * -----
    * x - initially contains v
    * g - temp data array for the algorithm
    * N - length of vectors x and g
    * matrix - length three array representing the tridiagonal stencil
    *
    * Output
    * ------
    * x - contains solution upon output (over-written)
    * g - is a working array, and will be modified as such
    **/

   int i;
   double m;

   g[0] = 0.0;  /* Assume the first row is the identity */

   /* loop from 1 to N - 2 inclusive, performing the forward sweep */
   for (i = 1; i < N - 1; i++)
   {
       m = 1.0 / (matrix[1] - matrix[0]*g[i-1]);
       g[i] = m*matrix[2];
       x[i] = m*(x[i] - matrix[0]*x[i-1]);
   }

   /* Do nothing for x[N-1], assume last row is the identity */

   /* loop from N - 2 to 1 inclusive to perform the back substitution */
   for (i = N - 2; i >= 1; i--)
   {
       x[i] = x[i] - g[i] * x[i+1];
   }
}

/* Helper function for Residual: Carry out mat-vec with tridiagonal stencil */
void
matvec_tridiag(PetscReal *x, PetscReal *g, PetscInt N, PetscReal* matrix)
{
   /**
    * Matvec solves g <-- Ax, where A is a tridiagonal matrix with stencil [ a, b, c].
    *
    * There is a built in assumption that the first and last rows are the
    * identity (boundary conditions)
    *
    * Input
    * -----
    * x - input vector
    * N - length of vectors x and g
    * matrix - length three array representing the tridiagonal stencil
    *
    * Output
    * ------
    * g - Equals A*x
    **/

   int i;

   /* loop from 1 to N - 2 inclusive, performing the matvec */
   for (i = 1; i < N - 1; i++)
   {
       g[i] = matrix[0]*x[i-1] + matrix[1]*x[i] + matrix[2]*x[i+1];
   }

   /* boundary points */
   g[0] = x[0];
   g[N-1] = x[N-1];

}

/* Compute three point backward Euler stencil */
void
compute_stencil(PetscReal   deltaX,
                PetscReal   deltaT,
                PetscReal * matrix)
{
   PetscReal cfl;

   cfl = (deltaT/(deltaX*deltaX));
   matrix[0] = -cfl;
   matrix[1] = 1.0 + 2*cfl;
   matrix[2] = -cfl;
}

/* Take a backward Euler step */
void
take_step(PetscReal * values,    /* state vector to evolve */
          PetscInt      size,      /* size of vector */
          PetscReal   t,         /* time value to evolve to */
          PetscReal   xstart,    /* left-most spatial grid point */
          PetscReal   deltaX,    /* spatial mesh size */
          PetscReal   deltaT,    /* temporal mesh size */
          PetscReal * matrix,    /* three-point matrix stencil */
          PetscReal * temp)      /* temporary array of size 'size' */
{
   int i;
   double x;
   /* Set up matrix stencil for 1D heat equation*/
   compute_stencil(deltaX, deltaT, matrix);

   /* Apply boundary conditions */
   values[0] = exact(t, 0.0);
   values[size-1] = exact(t, PI);

   /* PDE forcing */
   x = xstart;
   for(i = 0; i < size; i++)
   {
      values[i] = values[i] + deltaT*forcing(t, x);
      x = x + deltaX;
   }

   /* Backward Euler step */
   solve_tridiag(values, temp, size, matrix);

   Vec V;
   PetscReal Y;
   VecCreateSeqWithArray(PETSC_COMM_SELF, 1, size, values,&V);
   VecAssemblyBegin(V);
   VecAssemblyEnd(V);
   VecNorm(V,NORM_2,&Y);
   //printf("%f", Y);
}



/*--------------------------------------------------------------------------
 * Routines for coarsening and interpolating in space
 *--------------------------------------------------------------------------*/

/* Bilinear interpolation from grid size 2^{k-1} + 1
 * to grid size 2^{k} + 1 */
void
interpolate_1D(double * cvalues,
               double * fvalues,
               int      csize,
               int      fsize)
{
   int i;

   for (i = 1; i < fsize-1; i++)
   {
      if(i%2 == 1)
         (fvalues)[i] = 0.5*cvalues[i/2] + 0.5*cvalues[(i+1)/2];
      else
         (fvalues)[i] = cvalues[i/2];
   }

   /* Boundary Conditions */
   (fvalues)[0] = cvalues[0];
   (fvalues)[fsize-1] = cvalues[csize-1];
}

/* Bilinear coarsening from grid size 2^{k1} + 1
 * to grid size 2^{k-1} + 1 */
void
coarsen_1D(double * cvalues,
           double * fvalues,
           int      csize,
           int      fsize)
{
   int i, fidx;

   for (i = 1; i < csize-1; i++)
   {
      fidx = 2*i;
      cvalues[i] = 0.5*fvalues[fidx] + 0.25*fvalues[fidx+1] + 0.25*fvalues[fidx-1];
   }

   /* Boundary Conditions */
   cvalues[0] = fvalues[0];
   cvalues[csize-1] = fvalues[fsize-1];
}


/*--------------------------------------------------------------------------
 * Helper routines for output
 *--------------------------------------------------------------------------*/

/* Print accumulated info on space-time grids visited during the simulation */
void
print_sc_info(double *sc_info,
              int     max_levels)
{
   int i;
   double dx, dt;

   printf("\n-----------------------------------------------------------------\n");
   printf("-----------------------------------------------------------------\n\n");
   printf( " Per level diagnostic information \n\n");
   printf("level       dx          dt        dt/dx^2\n");
   printf("-----------------------------------------------------------------\n");
   for( i = 0; i < max_levels; i++)
   {
      dx = sc_info[i*2];
      dt = sc_info[i*2+1];
      if (dx == -1){
         break;
      }
      printf(" %2d   |   %1.2e    %1.2e    %1.2e\n", i, dx, dt, dt/(dx*dx) );
   }
   printf( "\n" );
}

/* Print a solution vector to file with the format
 *
 *    ntime_steps
 *    tstart
 *    tstop
 *    nspace_points
 *    xstart
 *    xstop
 *    x[0]
 *    x[1]
 *      .
 *      .
 *      .
 *    x[k]
 *
 **/
void
save_solution(char   *filename,
              double *values,
              int     size,
              double  xstart,
              double  xstop,
              int     ntime,
              double  tstart,
              double  tstop)
{
   FILE      *file;
   int i;

   file = fopen(filename, "w");
   fprintf(file, "%d\n",    ntime +1);
   fprintf(file, "%.14e\n", tstart );
   fprintf(file, "%.14e\n", tstop );
   fprintf(file, "%d\n",    size );
   fprintf(file, "%.14e\n", xstart );
   fprintf(file, "%.14e\n", xstop );

   for (i = 0; i < size; i++)
   {
      fprintf(file, "%.14e\n", values[i]);
   }

   fflush(file);
   fclose(file);
}

// /*--------------------------------------------------------------------------
//  * My App and Vector structures
//  *--------------------------------------------------------------------------*/
// /* can put anything in my app and name it anything as well */
// typedef struct _braid_App_struct
// {
//    MPI_Comm  comm;
//    double    tstart;       /* Define the temporal domain */
//    double    tstop;
//    int       ntime;
//    double    xstart;       /* Define the spatial domain */
//    double    xstop;
//    int       nspace;
//    double    matrix[3];    /* the three point spatial discretization stencil */
//    double *  g;            /* temporary vector for inversions and mat-vecs */
//    double *  sc_info;      /* Runtime information that tracks the space-time grids visited */
//    int       print_level;  /* Level of output desired by user (see the -help message below) */
// } my_App;
//
// /* Can put anything in my vector and name it anything as well */
// typedef struct _braid_Vector_struct
// {
//    PetscInt     size;
//    PetscReal *values;
//
// } my_Vector;
//
// /* create and allocate a vector */
// void
// create_vector(my_Vector **u,
//               int size)
// {
//    (*u) = (my_Vector *) malloc(sizeof(my_Vector));
//    ((*u)->size)   = size;
//    ((*u)->values) = (double *) malloc(size*sizeof(double));
// }
//
//
// /*--------------------------------------------------------------------------
//  * My integration routines
//  *--------------------------------------------------------------------------*/
//
// int my_Step(braid_App        app,
//             braid_Vector     ustop,
//             braid_Vector     fstop,
//             braid_Vector     u,
//             braid_StepStatus status)
// {
//    PetscReal tstart;             /* current time */
//    PetscReal tstop;              /* evolve to this time*/
//    PetscInt level, i;
//    PetscReal deltaX, deltaT;
//
//    braid_StepStatusGetLevel(status, &level);
//    braid_StepStatusGetTstartTstop(status, &tstart, &tstop);
//    deltaT = tstop - tstart;
//    deltaX = (app->xstop - app->xstart) / (ustop->size - 1.0);
//
//    /* XBraid forcing */
//    if(fstop != NULL)
//    {
//       for(i = 0; i < u->size; i++)
//       {
//          u->values[i] = u->values[i] + fstop->values[i];
//       }
//    }
//
//    /* Take backward Euler step
//     * Note: if an iterative solver were used, ustop->values would
//     *       contain the XBraid's best initial guess. */
//    take_step(u->values, u->size, tstop, app->xstart, deltaX, deltaT,
//          app->matrix, app->g);
//
//    /* Store info on space-time grids visited during the simulation */
//    (app->sc_info)[ (2*level) ] = deltaX;
//    (app->sc_info)[ (2*level) + 1] = deltaT;
//
//    /* no refinement */
//    braid_StepStatusSetRFactor(status, 1);
//
//    return 0;
// }
//
//
// int
// my_Init(braid_App     app,
//         double        t,
//         braid_Vector *u_ptr)
// {
//    my_Vector *u;
//    int    i;
//    int nspace = (app->nspace);
//    double deltaX = (app->xstop - app->xstart) / (nspace - 1.0);
//
//    /* Allocate vector */
//    create_vector(&u, nspace);
//
//    /* Initialize vector (with correct boundary conditions) */
//    if(t == 0.0)
//    {
//       /* Get the solution at time t=0 */
//       get_solution(u->values, u->size, 0.0, app->xstart, deltaX);
//    }
//    else
//    {
//       /* Use random values for u(t>0), this measures asymptotic convergence rate */
//       for(i=0; i < nspace; i++)
//       {
//          (u->values)[i] = ((double)braid_Rand())/braid_RAND_MAX;
//       }
//    }
//
//    *u_ptr = u;
//
//    return 0;
// }
//
// int
// my_Clone(braid_App     app,
//          braid_Vector  u,
//          braid_Vector *v_ptr)
// {
//    my_Vector *v;
//    int size = (u->size);
//    int i;
//
//    create_vector(&v, size);
//    for (i = 0; i < size; i++)
//    {
//       (v->values)[i] = (u->values)[i];
//    }
//    *v_ptr = v;
//
//    return 0;
// }
//
// int
// my_Free(braid_App    app,
//         braid_Vector u)
// {
//    free(u->values);
//    free(u);
//
//    return 0;
// }
//
// int
// my_Sum(braid_App     app,
//        double        alpha,
//        braid_Vector  x,
//        double        beta,
//        braid_Vector  y)
// {
//    int i;
//    int size = (y->size);
//
//    for (i = 0; i < size; i++)
//    {
//       (y->values)[i] = alpha*(x->values)[i] + beta*(y->values)[i];
//    }
//
//    return 0;
// }
//
// // int
// // my_Sum(braid_App     app,
// //        PetscReal        alpha,
// //        braid_Vector  x,
// //        PetscReal        beta,
// //        braid_Vector  y)
// // {
// //   Vec X,Y;
// //   VecCreateSeqWithArray(PETSC_COMM_SELF,x->size,x->size,x->values,&X);
// //   VecAssemblyBegin(X);
// //   VecAssemblyEnd(X);
// //   VecCreateSeqWithArray(PETSC_COMM_SELF,y->size,y->size,y->values,&Y);
// //   VecAssemblyBegin(Y);
// //   VecAssemblyEnd(Y);
// //   VecAXPBY(Y,alpha,beta,X);
// //
// //    return 0;
// // }
//
// int
// my_SpatialNorm(braid_App     app,
//                braid_Vector  u,
//                double       *norm_ptr)
// {
//    int    i;
//    int size   = (u->size);
//    double dot = 0.0;
//
//    for (i = 0; i < size; i++)
//    {
//       dot += (u->values)[i]*(u->values)[i];
//    }
//    *norm_ptr = sqrt(dot);
//
//    return 0;
// }
//
// int
// my_Access(braid_App          app,
//           braid_Vector       u,
//           braid_AccessStatus astatus)
// {
//    int        index, rank, level, done;
//    char       filename[255];
//    double     t, error;
//
//    braid_AccessStatusGetT(astatus, &t);
//    braid_AccessStatusGetTIndex(astatus, &index);
//    braid_AccessStatusGetLevel(astatus, &level);
//    braid_AccessStatusGetDone(astatus, &done);
//
//    /* Print solution to file if simulation is over */
//    if(done)
//    {
//       MPI_Comm_rank( (app->comm), &rank);
//       sprintf(filename, "%s.%07d.%05d", "ex-02.out", index, rank);
//       save_solution(filename, u->values, u->size, app->xstart,
//             app->xstop, app->ntime, app->tstart, app->tstop);
//    }
//
//    /* IF on the finest level AND print_level is high enough AND at the final time,
//     * THEN print out the discretization error */
//    if( (level == 0) && ((app->print_level) > 0) && (index == app->ntime) )
//    {
//       error = compute_error_norm(u->values, app->xstart, app->xstop, u->size, t);
//       printf("  Discretization error at final time:  %1.4e\n", error);
//       fflush(stdout);
//    }
//
//    return 0;
// }
//
// int
// my_BufSize(braid_App           app,
//            int                 *size_ptr,
//            braid_BufferStatus  bstatus)
// {
//    int size = (app->nspace);
//    *size_ptr = (size+1)*sizeof(double);
//    return 0;
// }
//
// int
// my_BufPack(braid_App           app,
//            braid_Vector        u,
//            void               *buffer,
//            braid_BufferStatus  bstatus)
// {
//    double *dbuffer = buffer;
//    int i, size = (u->size);
//
//    dbuffer[0] = size;
//    for (i = 0; i < size; i++)
//    {
//       dbuffer[i+1] = (u->values)[i];
//    }
//
//    braid_BufferStatusSetSize( bstatus,  (size+1)*sizeof(double));
//
//    return 0;
// }
//
// // int MyBraidApp::BufPack(braid_Vector       u_,
// //                           void               *buffer,
// //                           BraidBufferStatus  &status)
// // {
// //    BraidVector *u = (BraidVector*) u_;
// //    double *dbuffer = (double *) buffer;
// //
// //    dbuffer[0] = (u->value);
// //    status.SetSize(sizeof(double));
// //
// //    return 0;
// // }
//
// int
// my_BufUnpack(braid_App           app,
//              void               *buffer,
//              braid_Vector       *u_ptr,
//              braid_BufferStatus  bstatus)
// {
//    my_Vector *u = NULL;
//    double    *dbuffer = buffer;
//    int        i, size;
//
//    size = dbuffer[0];
//    create_vector(&u, size);
//
//    for (i = 0; i < size; i++)
//    {
//       (u->values)[i] = dbuffer[i+1];
//    }
//    *u_ptr = u;
//
//    return 0;
// }
//
//
// /*--------------------------------------------------------------------------
//  * my_Residual, my_Coarsen and my_Refine are advanced XBraid options, ignore
//  * them until you understand the rest of the driver.
//  *--------------------------------------------------------------------------*/
//
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
//
//    braid_StepStatusGetTstartTstop(status, &tstart, &tstop);
//    deltaT = tstop - tstart;
//    deltaX = (app->xstop - app->xstart) / (ustop->size - 1.0);
//
//    /* Set up matrix stencil for 1D heat equation*/
//    compute_stencil(deltaX, deltaT, app->matrix);
//
//    /* Residual r = A*xstop - r - forcing - boundary
//     *   note: there are no boundary terms here */
//    matvec_tridiag(ustop->values, app->g, ustop->size, app->matrix);
//    x = app->xstart;
//    for(i = 0; i < r->size; i++)
//    {
//       r->values[i] = app->g[i] - r->values[i] - deltaT*forcing(tstop, x);
//       x = x + deltaX;
//    }
//
//    return 0;
// }
//
// /* Bilinear Coarsening */
// int
// my_Coarsen(braid_App              app,
//            braid_Vector           fu,
//            braid_Vector          *cu_ptr,
//            braid_CoarsenRefStatus status)
// {
//
//    int csize, level;
//    my_Vector *v;
//
//    /* This returns the level for fu */
//    braid_CoarsenRefStatusGetLevel(status, &level);
//
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
//
//    *cu_ptr = v;
//
//    return 0;
// }
//
// /* Bilinear interpolation */
// int
// my_Interp(braid_App              app,
//           braid_Vector           cu,
//           braid_Vector          *fu_ptr,
//           braid_CoarsenRefStatus status)
// {
//
//    int fsize, level;
//    my_Vector *v;
//
//    /* This returns the level for fu_ptr */
//    braid_CoarsenRefStatusGetLevel(status, &level);
//
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
//
//    *fu_ptr = v;
//
//    return 0;
// }
