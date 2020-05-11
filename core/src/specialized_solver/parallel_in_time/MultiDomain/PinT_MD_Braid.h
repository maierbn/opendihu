// specific functions for XBraid with Implicit Euler
// defines my_Step (integration step)

#pragma once

#include <braid.h>

//#include "time_stepping_scheme/implicit_euler.h"
#include "time_stepping_scheme/heun.h"
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

int my_Step_MD(braid_App        app,
            braid_Vector     ustop,
            braid_Vector     fstop,
            braid_Vector     u,
            braid_StepStatus status);

int
my_Init_MD(braid_App     app,
        double        t,
        braid_Vector *u_ptr);