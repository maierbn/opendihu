#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

// #include "braid.hpp"
#include "PinT_IE_lib.c"

// --------------------------------------------------------------------------
// User-defined routines and objects
// --------------------------------------------------------------------------

// Define BraidVector, can contain anything, and be named anything
// --> Put all time-dependent information here
class BraidVector
{
public:
   // Each vector holds the scalar solution value at a particular time
   double value;

   // Construct a BraidVector for a given double
   BraidVector(double value_) : value(value_) { }

   // Deconstructor
   virtual ~BraidVector() {};

};


// Wrapper for BRAID's App object
// --> Put all time INDEPENDENT information here
class MyBraidApp : public BraidApp
{
protected:
   // BraidApp defines tstart, tstop, ntime and comm_t

public:
   // Constructor
   MyBraidApp(MPI_Comm comm_t_, int rank_, double tstart_ = 0.0, double tstop_ = 1.0, int ntime_ = 100);

   // We will need the MPI Rank
   int rank;

   // Deconstructor
   virtual ~MyBraidApp() {};

   // Define all the Braid Wrapper routines
   // Note: braid_Vector == BraidVector*
   virtual int Step(braid_Vector    u_,
                    braid_Vector    ustop_,
                    braid_Vector    fstop_,
                    BraidStepStatus &pstatus);

   virtual int Clone(braid_Vector  u_,
                     braid_Vector *v_ptr);

   virtual int Init(double        t,
                    braid_Vector *u_ptr);

   virtual int Free(braid_Vector u_);

   virtual int Sum(double       alpha,
                   braid_Vector x_,
                   double       beta,
                   braid_Vector y_);

   virtual int SpatialNorm(braid_Vector  u_,
                           double       *norm_ptr);

   virtual int BufSize(int *size_ptr,
                       BraidBufferStatus  &status);

   virtual int BufPack(braid_Vector  u_,
                       void         *buffer,
                       BraidBufferStatus  &status);

   virtual int BufUnpack(void         *buffer,
                         braid_Vector *u_ptr,
                         BraidBufferStatus  &status);

   virtual int Access(braid_Vector       u_,
                      BraidAccessStatus &astatus);

   // Not needed in this example
   virtual int Residual(braid_Vector     u_,
                        braid_Vector     r_,
                        BraidStepStatus &pstatus) { return 0; }

   // Not needed in this example
   virtual int Coarsen(braid_Vector   fu_,
                       braid_Vector  *cu_ptr,
                       BraidCoarsenRefStatus &status) { return 0; }

   // Not needed in this example
   virtual int Refine(braid_Vector   cu_,
                      braid_Vector  *fu_ptr,
                      BraidCoarsenRefStatus &status)  { return 0; }

};

// Braid App Constructor
MyBraidApp::MyBraidApp(MPI_Comm comm_t_, int rank_, double tstart_, double tstop_, int ntime_)
   : BraidApp(comm_t_, tstart_, tstop_, ntime_)
{
   rank = rank_;
}

//
int MyBraidApp::Step(braid_Vector    u_,
                     braid_Vector    ustop_,
                     braid_Vector    fstop_,
                     BraidStepStatus &pstatus)
{

   BraidVector *u = (BraidVector*) u_;
   double tstart;             // current time
   double tstop;              // evolve to this time

   // Get time step information
   pstatus.GetTstartTstop(&tstart, &tstop);

   // Use backward Euler to propagate solution
   (u->value) = 1./(1. + tstop-tstart)*(u->value);

   // no refinement
   pstatus.SetRFactor(1);

   return 0;
   // BraidVector *u = (BraidVector*) u_;
   // PetscReal tstart;             /* current time */
   // PetscReal tstop;              /* evolve to this time*/
   // PetscInt level, i;
   // PetscReal deltaX, deltaT;
   //
   // pstatus.GetLevel(&level);
   // pstatus.GetTstartTstop(&tstart, &tstop);
   // deltaT = tstop - tstart;
   // deltaX = (app->xstop - app->xstart) / (ustop->size - 1.0);
   //
   // /* XBraid forcing */
   // if(fstop != NULL)
   // {
   //    for(i = 0; i < u->size; i++)
   //    {
   //       u->values[i] = u->values[i] + fstop->values[i];
   //    }
   // }
   //
   // /* Take backward Euler step
   //  * Note: if an iterative solver were used, ustop->values would
   //  *       contain the XBraid's best initial guess. */
   // take_step(u->values, u->size, tstop, app->xstart, deltaX, deltaT,
   //       app->matrix, app->g);
   //
   // /* Store info on space-time grids visited during the simulation */
   // (app->sc_info)[ (2*level) ] = deltaX;
   // (app->sc_info)[ (2*level) + 1] = deltaT;
   //
   // /* no refinement */
   // braid_StepStatusSetRFactor(status, 1);
   //
   // return 0;

}

int MyBraidApp::Init(double        t,
                       braid_Vector *u_ptr)
{
   BraidVector *u = new BraidVector(0.0);
   if (t != tstart)
   {
      u->value = 0.456;
   }
   else
   {
      u->value = 1.0;
   }

   *u_ptr = (braid_Vector) u;
   return 0;

}

int MyBraidApp::Clone(braid_Vector  u_,
                        braid_Vector *v_ptr)
{
   BraidVector *u = (BraidVector*) u_;
   BraidVector *v = new BraidVector(u->value);
   *v_ptr = (braid_Vector) v;

   return 0;
}


int MyBraidApp::Free(braid_Vector u_)
{
   BraidVector *u = (BraidVector*) u_;
   delete u;
   return 0;
}

int MyBraidApp::Sum(double       alpha,
                      braid_Vector x_,
                      double       beta,
                      braid_Vector y_)
{
   BraidVector *x = (BraidVector*) x_;
   BraidVector *y = (BraidVector*) y_;
   (y->value) = alpha*(x->value) + beta*(y->value);
   return 0;
}

int MyBraidApp::SpatialNorm(braid_Vector  u_,
                              double       *norm_ptr)
{
   double dot;
   BraidVector *u = (BraidVector*) u_;
   dot = (u->value)*(u->value);
   *norm_ptr = sqrt(dot);
   return 0;
}

int MyBraidApp::BufSize(int                *size_ptr,
                          BraidBufferStatus  &status)
{
   *size_ptr = sizeof(double);
   return 0;
}

int MyBraidApp::BufPack(braid_Vector       u_,
                          void               *buffer,
                          BraidBufferStatus  &status)
{
   BraidVector *u = (BraidVector*) u_;
   double *dbuffer = (double *) buffer;

   dbuffer[0] = (u->value);
   status.SetSize(sizeof(double));

   return 0;
}

int MyBraidApp::BufUnpack(void              *buffer,
                            braid_Vector      *u_ptr,
                            BraidBufferStatus &status)
{
   double *dbuffer = (double *) buffer;

   BraidVector *u = new BraidVector(dbuffer[0]);
   *u_ptr = (braid_Vector) u;

   return 0;
}

int MyBraidApp::Access(braid_Vector       u_,
                         BraidAccessStatus &astatus)
{
   char       filename[255];
   FILE      *file;
   BraidVector *u = (BraidVector*) u_;

   // Extract information from astatus
   int done, level, iter, index;
   double t;
   astatus.GetTILD(&t, &iter, &level, &done);
   // astatus.GetTIndex(&index);

   // Print information to file
   sprintf(filename, "%s.%04d.%03d", "ex-01.out", index, rank);
   file = fopen(filename, "w");
   fprintf(file, "%.14e\n", (u->value));
   fflush(file);
   fclose(file);

   return 0;
}
