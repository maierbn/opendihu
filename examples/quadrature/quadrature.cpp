#include <iostream>
#include <cstdlib>
#include <cmath>

#include "opendihu.h"

// declaration of the test function which is defined in nonlinear_test.cpp
void testQuadrature(int argc, char *argv[]);

// function to be integrated numerically
double function1(double x)
{
  return cos(x)*x;
}

void testGauss2()
{
  Quadrature::Gauss<2> gauss2;
  
  // get the number of evaluations which is the number of sampling points, i.e. 2 in this case
  const int nEvaluations = gauss2.numberEvaluations();
  std::cout << "number evaluations: " << nEvaluations << std::endl;
  
  // get sampling points, ie "Stützstellen"
  std::array<double,nEvaluations> samplingPoints = gauss2.samplingPoints();
  
  // evaluate function at sampling points
  std::array<double,nEvaluations> evaluations({0.0});
  
  // loop over sampling points
  for (int i = 0; i < nEvaluations; i++)
  {
    double samplingPoint = samplingPoints[i];
    evaluations[i] = function1(samplingPoint);
  }
  
  // perform the quadrature in [0,1]
  double result = gauss2.computeIntegral(evaluations);
  
  // compute reference solution and error
  double referenceSolution = -1 + sin(1) + cos(1);    //wolfram alpha: integrate(cos(x)*x,x,0,1)
  double error = fabs(result - referenceSolution);
  
  // output result
  std::cout << "Result: " << result 
    << ", error: " << error << std::endl;
}

int main(int argc, char *argv[])
{
  testGauss2();
  
  testQuadrature(argc, argv);
  
  return EXIT_SUCCESS;
}