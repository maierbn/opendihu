#pragma once

#include <functional>

// source: https://github.com/matteotiziano/nelder-mead
// author: Matteo, user matteotiziano
// licence: MIT

namespace MathUtility
{

namespace NelderMead
{

//-----------------------------------------------------------------------------
// Definitions
//-----------------------------------------------------------------------------

// define a generic point containing a position (x) and a value (fx)
typedef struct
{
  double *x;
  double fx;
} point_t;

// define a simplex struct containing an array of n+1 points (p)
// each having dimension (n)
typedef struct
{
  point_t *p;
  int n;
} simplex_t;

// define optimization settings
typedef struct
{
  double tolx;
  double tolf;
  int max_iter;
  int max_eval;
  int verbose;
} optimset_t;

//-----------------------------------------------------------------------------
// Cost function interface
//-----------------------------------------------------------------------------

typedef std::function<void (int, point_t *, const void *)> fun_t;


//-----------------------------------------------------------------------------
// Main function
// - n is the dimension of the data
// - start is the initial point (unchanged in output)
// - solution is the minimizer
// - cost_function is a pointer to a fun_t type function
// - args are the optional arguments of cost_function
// - optimset are the optimisation settings
//-----------------------------------------------------------------------------
void optimize(int n, const point_t *start, point_t *solution,
              fun_t cost_function, const void *args,
              const optimset_t *optimset);

//-----------------------------------------------------------------------------
// Utility functions
//-----------------------------------------------------------------------------

int compare(const void *, const void *);

void simplexSort(simplex_t *);

void getCentroid(const simplex_t *, point_t *);

double modulus(double);

int continueMinimization(const simplex_t *, int, int, const optimset_t *);

void updatePoint(const simplex_t *, const point_t *, double, point_t *);

void copyPoint(int, const point_t *, point_t *);

void swapPoints(int, point_t *, point_t *);

}  // namespace NelderMead
}  // namespace MathUtility
