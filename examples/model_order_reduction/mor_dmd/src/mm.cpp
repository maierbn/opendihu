#include <iostream>
#include <cstdlib>
#include <vector>
#include "opendihu.h"
#include "utility/svd_utility.h"
#include "utility/dmd_utility.h"
#include <complex>

int main(int argc, char *argv[])
{
	/*
	double six[] {
		2, 3, 5, 7, 11, 13
	};

	double four[] {
		2, 3, 5, 7
	};

	double* fourTransposed = new double[4];
*/
/*double ten[] {
	2, 3, 5, 7, 11, 13,
	17, 19, 23, 29
};*/

/*double fifteen[] {
	2, 3, 5, 7, 11, 13,
	17, 19, 23, 29,
	31, 37, 41, 43, 47
};*/

//double* i = new double[10];
//double ii = new double[15];
//double iii = new double[6];
//double iv = new double[6];
//double v = new double[15];
//double vi = new double[10];
/*
	double* test = new double[6];

	SvdUtility::printMatrix("four", four, 2, 2);

	std::cout << "four[2] = " << four[2] << endl;

	SvdUtility::getMatrixProduct(six, four, test, 3, 2, 2, false);

	SvdUtility::printMatrix("four", four, 2, 2);

	SvdUtility::printMatrix("test", test, 3, 2);

	SvdUtility::printMatrix("four", four, 2, 2);

	SvdUtility::transposeMatrix(four, fourTransposed, 2, 2);

	SvdUtility::printMatrix("four", fourTransposed, 2, 2);

	SvdUtility::getMatrixProduct(six, fourTransposed, test, 3, 2, 2, false);

	SvdUtility::printMatrix("test", test, 3, 2);


	int j0 = 3;
	int k0 = 4;
	double v0[] {
		1, -1,  0,
		0, -2,  1,
		1,  0, -1,
		2,  0,  0
	};

	// n=length(sigmas);
	int n0 = std::min(j0, k0);

	// J > K => only the first K cols of U are computed
	double* u0 = new double[j0 * n0];

	// J < K => only the first J rows of T^T are computed
	double* tTransposed0 = new double[n0 * k0];

	double* sigmas0 = new double[n0];
	double* sigma0 = new double[n0 * n0];

	// [U,Sigma,T]=svd(V,'econ');
	SvdUtility::getSVD(v0, j0, k0, u0, sigmas0, tTransposed0, sigma0);

	// print resulting matrices
	SvdUtility::printMatrix("U0", u0, j0, n0);
	SvdUtility::printMatrix("Sigma0", sigma0, n0, n0);
	SvdUtility::printMatrix("T0^T", tTransposed0, n0, k0);

	int j1 = 4;
	int k1 = 3;
	double v1[] {
	1, -1,  0, 2,
	0, -2,  1, 0,
	1,  0, -1, 0
	};

	// n=length(sigmas);
	int n1 = std::min(j1, k1);

	// J > K => only the first K cols of U are computed
	double* u1 = new double[j1 * n1];

	// J < K => only the first J rows of T^T are computed
	double* tTransposed1 = new double[n1 * k1];

	double* sigmas1 = new double[n1];
	double* sigma1 = new double[n1 * n1];

	// [U,Sigma,T]=svd(V,'econ');
	SvdUtility::getSVD(v1, j1, k1, u1, sigmas1, tTransposed1, sigma1);

	// print resulting matrices
	SvdUtility::printMatrix("U1", u1, j1, n1);
	SvdUtility::printMatrix("Sigma", sigma1, n1, n1);
	SvdUtility::printMatrix("T'", tTransposed1, n1, k1);
	*/

	/*double _Complex* matrixA = new double _Complex[3 * 3];
	for (int col = 0; col < 3; ++col)
	{
		for (int row = 0; row < 3; ++row)
		{
			matrixA[row + col * 3] = csqrt(-1.0) * row + col;
		}
	}*/

	/*double _Complex* matrixB = new double _Complex[3 * 5];
	for (int col = 0; col < 5; ++col)
	{
		for (int row = 0; row < 3; ++row)
		{
			matrixB[row + col * 3] = csqrt(-1.0) * col + row;
		}
	}*/
	// SvdUtility::printMatrix("matrixA", matrixA, 3, 3);
	//SvdUtility::printMatrix("matrixB", matrixB, 3, 5);

	/*double _Complex* matrixC = new double _Complex[3 * 3];
	for (int col = 0; col < 3; ++col)
	{
		for (int row = 0; row < 3; ++row)
		{
			matrixC[row + col * 3] = 0;
		}
	}
	SvdUtility::printMatrix("matrixC", matrixC, 3, 3);
	SvdUtility::getMatrixPower(matrixA, matrixC, 3, 0);
	SvdUtility::printMatrix("matrixC", matrixC, 3, 3);*/

	/*double _Complex negOne = csqrt(-1.0) * csqrt(-1.0);
	std::cout << creal(negOne) << " + " << cimag(negOne) << "i" << endl;*/

	int j = 4;
	//int k = 4;
	double v[]{
		1, -1,  0, 2,
		0, -2,  1, 0,
		1,  0, -1, 0,
		0,  1,  0, 0
	};

	double _Complex* ev = new double _Complex[j];
	double _Complex* lev = new double _Complex[j * j];

	DmdUtility::getEigen(v, j, ev, lev);

	DmdUtility::printMatrix("ev", ev, 1, j);
	DmdUtility::printMatrix("lev", lev, j, j);

	/*for (int i = 0; i < j; ++i)
	{
		for (int k = 0; k < j; ++k)
		{
			std::cout << lev[2 * k + 2 * i * j] << " + " << lev[2 * k + 2 * i * j + 1] << "i ";
		}
		std::cout << endl;
	}*/

	return EXIT_SUCCESS;
}
