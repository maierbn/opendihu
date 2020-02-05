#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

// int MyBraidApp::BufPack(braid_Vector       u_,
//                           void               *buffer,
//                           BraidBufferStatus  &status)
int my_BufPack(braid_Vector       u_,
                          void               *buffer,
                          BraidBufferStatus  &status)
{
   BraidVector *u = (BraidVector*) u_;
   double *dbuffer = (double *) buffer;

   dbuffer[0] = (u->value);
   status.SetSize(sizeof(double));

   return 0;
}

// int MyBraidApp::BufUnpack(void              *buffer,
//                             braid_Vector      *u_ptr,
//                             BraidBufferStatus &status)
int my_BufUnpack(void              *buffer,
                            braid_Vector      *u_ptr,
                            BraidBufferStatus &status)
{
   double *dbuffer = (double *) buffer;

   BraidVector *u = new BraidVector(dbuffer[0]);
   *u_ptr = (braid_Vector) u;

   return 0;
}
