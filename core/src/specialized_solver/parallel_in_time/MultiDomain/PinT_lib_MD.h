#pragma once

#include <braid.h>

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include <petscdm.h>
#include <petscdmda.h>
#include <petscts.h>
#include <petscdraw.h>

//#include "specialized_solver/parallel_in_time/ImplicitEuler/PinT_IE_Braid.h"
//#include "specialized_solver/parallel_in_time/MultiDomain/PinT_MD_Braid.h"

#ifdef M_PI
   #define PI M_PI
#else
   #define PI 3.14159265358979
#endif

/*--------------------------------------------------------------------------
 * Routines defining solution, forcing term and how to compute the error
 *--------------------------------------------------------------------------*/

/* Exact solution */
double exact(PetscReal t, PetscReal x);

/* Initialize array of values to solution at time t */
void
get_solution(PetscReal * values,
             PetscInt      size,
             PetscReal   t,
             PetscReal   xstart,
             PetscReal   deltaX);

/* Forcing term F(t, x) for PDE,  u_t = u_xx + F(t,x) */
double forcing(PetscReal t, PetscReal x);

/* Compute L2-norm of the error at a point in time */
double compute_error_norm(double * values,
                          double   xstart,
                          double   xstop,
                          int      nspace,
                          double   t);


/*--------------------------------------------------------------------------
 * Routines for computing time-steps and residuals
 *--------------------------------------------------------------------------*/

/* Helper function for tridiagonal solver */
PetscReal dabs(PetscReal x);


/* Helper function for Step: Tridiagonal system solver (Thomas algorithm) */
void
solve_tridiag(PetscReal *x, PetscReal *g, PetscInt N, PetscReal* matrix);

/* Helper function for Residual: Carry out mat-vec with tridiagonal stencil */
void
matvec_tridiag(PetscReal *x, PetscReal *g, PetscInt N, PetscReal* matrix);

/* Compute three point backward Euler stencil */
void
compute_stencil(PetscReal   deltaX,
                PetscReal   deltaT,
                PetscReal * matrix);

/* Take a backward Euler step */
void
take_step(PetscReal * values,    /* state vector to evolve */
          PetscInt      size,      /* size of vector */
          PetscReal   t,         /* time value to evolve to */
          PetscReal   xstart,    /* left-most spatial grid point */
          PetscReal   deltaX,    /* spatial mesh size */
          PetscReal   deltaT,    /* temporal mesh size */
          PetscReal * matrix,    /* three-point matrix stencil */
          PetscReal * temp);      /* temporary array of size 'size' */


/*--------------------------------------------------------------------------
 * Routines for coarsening and interpolating in space
 *--------------------------------------------------------------------------*/

/* Bilinear interpolation from grid size 2^{k-1} + 1
 * to grid size 2^{k} + 1 */
void
interpolate_1D(double * cvalues,
               double * fvalues,
               int      csize,
               int      fsize);

/* Bilinear coarsening from grid size 2^{k1} + 1
 * to grid size 2^{k-1} + 1 */
void
coarsen_1D(double * cvalues,
           double * fvalues,
           int      csize,
           int      fsize);


/*--------------------------------------------------------------------------
 * Helper routines for output
 *--------------------------------------------------------------------------*/

/* Print accumulated info on space-time grids visited during the simulation */
void
print_sc_info(double *sc_info,
              int     max_levels);

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
              double  tstop);