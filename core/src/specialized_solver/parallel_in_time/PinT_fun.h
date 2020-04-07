// most important parallel-in-time structures and functions for XBraid
// more functions in PinT_lib.h

#pragma once

#include <braid.h>

#include "time_stepping_scheme/implicit_euler.h"
#include "spatial_discretization/finite_element_method/finite_element_method.h"
#include "basis_function/lagrange.h"
#include "mesh/structured_regular_fixed.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>

#include <petscdm.h>
#include <petscdmda.h>
#include <petscts.h>
#include <petscdraw.h>
#include <petscvec.h>

/*--------------------------------------------------------------------------
 * My App and Vector structures
 *--------------------------------------------------------------------------*/
/* can put anything in my app and name it anything as well */
typedef struct _braid_App_struct
{
  MPI_Comm  comm;
  double    tstart;       /* Define the temporal domain */
  double    tstop;
  int       ntime;
  double    xstart;       /* Define the spatial domain */
  double    xstop;
  int       nspace;
  double    matrix[3];    /* the three point spatial discretization stencil */
  double *  g;            /* temporary vector for inversions and mat-vecs */
  double *  sc_info;      /* Runtime information that tracks the space-time grids visited */
  int       print_level;  /* Level of output desired by user (see the -help message below) */

  typedef  TimeSteppingScheme::ImplicitEuler<
    SpatialDiscretization::FiniteElementMethod<
      Mesh::StructuredRegularFixedOfDimension<1>,
      BasisFunction::LagrangeOfOrder<>,
      Quadrature::None,
      Equation::Dynamic::IsotropicDiffusion
    >
  > NestedSolver;

  std::vector<
    std::shared_ptr<
      NestedSolver
    >
  > *implicitEulerSolvers;   //< vector of nested solvers (implicit euler) for solution on different grids

} my_App;

/* Can put anything in my vector and name it anything as well */
typedef struct _braid_Vector_struct
{
   PetscInt     size;
   PetscReal *values;

} my_Vector;

/* create and allocate a vector */
void
create_vector(my_Vector **u,
              int size);

int
my_Init(braid_App     app,
        double        t,
        braid_Vector *u_ptr);

int
my_Clone(braid_App     app,
         braid_Vector  u,
         braid_Vector *v_ptr);

int
my_Free(braid_App    app,
        braid_Vector u);

int
my_Sum(braid_App     app,
       double        alpha,
       braid_Vector  x,
       double        beta,
       braid_Vector  y);

int
my_SpatialNorm(braid_App     app,
               braid_Vector  u,
               double       *norm_ptr);

int
my_Access(braid_App          app,
          braid_Vector       u,
          braid_AccessStatus astatus);

int
my_BufSize(braid_App           app,
           int                 *size_ptr,
           braid_BufferStatus  bstatus);

int
my_BufPack(braid_App           app,
           braid_Vector        u,
           void               *buffer,
           braid_BufferStatus  bstatus);

int
my_BufUnpack(braid_App           app,
             void               *buffer,
             braid_Vector       *u_ptr,
             braid_BufferStatus  bstatus);

int
my_Residual(braid_App        app,
            braid_Vector     ustop,
            braid_Vector     r,
            braid_StepStatus status);

/* Bilinear Coarsening */
int
my_Coarsen(braid_App              app,
           braid_Vector           fu,
           braid_Vector          *cu_ptr,
           braid_CoarsenRefStatus status);

/* Bilinear interpolation */
int
my_Interp(braid_App              app,
          braid_Vector           cu,
          braid_Vector          *fu_ptr,
          braid_CoarsenRefStatus status);
