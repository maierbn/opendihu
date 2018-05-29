#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <limits>

#include "opendihu.h"

#define numNodes 5
#define _PI (2.0 * asin(1.0))

// declaration of the test function which is defined in nonlinear_test.cpp
// void testQuadrature(int argc, char *argv[]);

// function to be integrated numerically
double function1(double x) {
  return cos(x)*x;
}

double function2(double x) {
	return (x*exp(x)) / ((x + 1.)*(x + 1.));
}

double function3(double x) {
	return 2.*x*x*x + 5.*x*x + 0.2*x + 9;
}

void testGauss()
{
  Quadrature::Gauss<numNodes> gauss;
  
  // get the number of evaluations which is the number of sampling points, i.e. 2 in this case
  const int nEvaluations = numNodes;// gauss2.numberEvaluations();
  std::cout << "number evaluations: " << nEvaluations << std::endl;
  
  // get sampling points, ie "Stützstellen"
  std::array<double,nEvaluations> samplingPoints = gauss.samplingPoints();
  
  // evaluate function at sampling points
  std::array<double,nEvaluations> evaluations{0.0};
  
  // loop over sampling points
  for (int i = 0; i < nEvaluations; i++)
  {
    double samplingPoint = samplingPoints[i];
    evaluations[i] = function2(samplingPoint);
  }
  
  // perform the quadrature in [0,1]
  double result = gauss.computeIntegral(evaluations);
  
  // compute reference solution and error
  //double referenceSolution = -1 + sin(1) + cos(1);    //wolfram alpha: integrate(cos(x)*x,x,0,1)
  //double error = fabs(result - referenceSolution);

  double referenceSolution = (exp(1)-2.)/2.;    //wolfram alpha: int (xe^x)/(x+1)^2 dx from 0 to 1
  double error = fabs(result - referenceSolution);

  //double referenceSolution = 11.2667;    //wolfram alpha: int 2.*x*x*x + 5.*x*x + 0.2*x + 9 dx from 0 to 1
  //double error = fabs(result - referenceSolution);
  
  // output result
  std::cout << "Result Gauss: " << result 
    << ", error: " << error << std::endl;
}

void testNC()
{
	Quadrature::NC<numNodes> nc;

	// get the number of evaluations which is the number of sampling points, i.e. 2 in this case
	const int nEvaluations = numNodes;// gauss2.numberEvaluations();
	std::cout << "number evaluations: " << nEvaluations << std::endl;

	// get sampling points, ie "Stützstellen"
	std::array<double, nEvaluations> samplingPoints = nc.samplingPoints();

	// evaluate function at sampling points
	std::array<double, nEvaluations> evaluations{ 0.0 };

	// loop over sampling points
	for (int i = 0; i < nEvaluations; i++)
	{
		double samplingPoint = samplingPoints[i];
		evaluations[i] = function2(samplingPoint);
	}

	// perform the quadrature in [0,1]
	double result = nc.computeIntegral(evaluations);

	// compute reference solution and error
	//double referenceSolution = -1 + sin(1) + cos(1);    //wolfram alpha: integrate(cos(x)*x,x,0,1)
	//double error = fabs(result - referenceSolution);

	double referenceSolution = (exp(1)-2.)/2.;    //wolfram alpha: int (xe^x)/(x+1)^2 dx from 0 to 1
	double error = fabs(result - referenceSolution);

	//double referenceSolution = 11.2667;    //wolfram alpha: int 2.*x*x*x + 5.*x*x + 0.2*x + 9 dx from 0 to 1
	//double error = fabs(result - referenceSolution);

	// output result
	std::cout << "Result Newton-Cotes: " << result
		<< ", error: " << error << std::endl;
}

void testCC()
{
	Quadrature::CC<numNodes> cc;

	// get the number of evaluations which is the number of sampling points, i.e. 2 in this case
	const int nEvaluations = numNodes;// gauss2.numberEvaluations();
	std::cout << "number evaluations: " << nEvaluations << std::endl;

	// get sampling points, ie "Stützstellen"
	std::array<double, nEvaluations> samplingPoints = cc.samplingPoints();

	// evaluate function at sampling points
	std::array<double, nEvaluations> evaluations{ 0.0 };

	// loop over sampling points
	for (int i = 0; i < nEvaluations; i++)
	{
		double samplingPoint = samplingPoints[i];
		evaluations[i] = function2(samplingPoint);
	}

	// perform the quadrature in [0,1]
	double result = cc.computeIntegral(evaluations);

	// compute reference solution and error
	//double referenceSolution = -1 + sin(1) + cos(1);    //wolfram alpha: integrate(cos(x)*x,x,0,1)
	//double error = fabs(result - referenceSolution);

	double referenceSolution = (exp(1)-2.)/2.;    //wolfram alpha: int (xe^x)/(x+1)^2 dx from 0 to 1
	double error = fabs(result - referenceSolution);

	//double referenceSolution = 11.2667;    //wolfram alpha: int 2.*x*x*x + 5.*x*x + 0.2*x + 9 dx from 0 to 1
	//double error = fabs(result - referenceSolution);

	// output result
	std::cout << "Result Clenshaw-Curtis: " << result
		<< ", error: " << error << std::endl;
}

void CCweights() {
	#define numNodesCC 9
	std::vector<double> weights;
	for (int i = 0; i < numNodesCC; ++i) {
		double coef{};
		for (int j = 1; j <= (numNodesCC - 1) / 2.0; ++j) {
			double frac = (1.0 / (4.0 * j * j - 1.0));
			double cosine = 2.0 * cos((2.0 * _PI * i * j) / static_cast<double>(numNodesCC - 1));
			if (2.0 * j == numNodesCC - 1) cosine /= 2.0;
			coef += frac * cosine;
		}
		double w = (2.0 / static_cast<double>(numNodesCC - 1)) * (1.0 - coef);
		weights.push_back(w);
	}
	std::cout.precision(16);
	for (int i = 0; i <= numNodesCC - 1; ++i) {
		if (i == 0 || i == numNodesCC - 1) weights[i] /= 2.;
		// std::cout << weights[i] << "," << std::endl;
	}
}

int main(int argc, char *argv[])
{
  testGauss();
  testNC();
  testCC();
  
  //testQuadrature(argc, argv);
  std::cin.get();
  return EXIT_SUCCESS;
}