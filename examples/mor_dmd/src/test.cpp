#include "opendihu.h"
#include "utility/svd_utility.h"

#include <iostream>
#include <cstdlib>
#include <vector>
#include <complex>
#include <math.h>

// using namespace std;

int main(int argc, char *argv[])
{
	// input data
	//std::vector<double> parsedCsv = SvdUtility::readCSV("./out/data.csv");
	//double v[parsedCsv.size()];
	//copy(parsedCsv.begin(), parsedCsv.end(), v);

	//int rowCount = SvdUtility::getCSVRowCount("./out/data.csv");
	//int columnCount = SvdUtility::getCSVColumnCount("./out/data.csv");

	double varepsilon1 = 0.1;
	double deltat = 0.01;

	// simple square matrix
	/*int j = 3;
	int k = 3;
	double v[] {
	1, -1,  0,
	0, -2,  1,
	1,  0, -1
	};*/

	// simple wide and short matrix
	/*int j = 3;
	int k = 4;
	double v[] {
		1, -1,  0,
		0, -2,  1,
		1,  0, -1,
		2,  0,  0
	};*/

	// simple high and narrow matrix
	/*int j = 4;
	int k = 3;
	double v[] {
	1, -1,  0, 2,
	0, -2,  1, 0,
	1,  0, -1, 0
	};*/

	// [J,K]=size(V);

	// complex high and narrow matrix
	/*int j = 6;
	int k = 5;
	double v[] {
			  8.79,  6.11, -9.15,  9.57, -3.49,  9.84,
			  9.93,  6.91, -7.93,  1.64,  4.02,  0.15,
			  9.83,  5.04,  4.86,  8.83,  9.80, -8.99,
			  5.45, -0.27,  4.85,  0.74, 10.00, -6.02,
			  3.16,  7.98,  3.01,  5.80,  4.27, -5.31
	};*/

	// complex square matrix
	/*int j = 5;
	int k = 5;
	double v[j * k] = {
			  8.79,  6.11, -9.15,  9.57, -3.49,
			  9.93,  6.91, -7.93,  1.64,  4.02,
			  9.83,  5.04,  4.86,  8.83,  9.80,
			  5.45, -0.27,  4.85,  0.74, 10.00,
			  3.16,  7.98,  3.01,  5.80,  4.27
	};*/

	// complex wide and short matrix
	int j = 5;
	int k = 6;
	double v[]{
			  8.79,  6.11, -9.15,  9.57, -3.49,
			  9.93,  6.91, -7.93,  1.64,  4.02,
			  9.83,  5.04,  4.86,  8.83,  9.80,
			  5.45, -0.27,  4.85,  0.74, 10.00,
			  3.16,  7.98,  3.01,  5.80,  4.27,
			  9.84,  0.15, -8.99, -6.02, -5.31
	};

	// print input matrix sizes, j rows, k cols
	std::cout << "J = rows = " << j << endl << "K = cols =  " << k << endl << endl;
	SvdUtility::printMatrix("V", v, j, k);

	// n=length(sigmas);
	int n = std::min(j, k);

	// J > K => only the first K cols of U are computed
	double* u = new double[j * n];

	// J < K => only the first J rows of T^T are computed
	double* tTransposed = new double[n * k];

	double* sigmas = new double[n];
	double* sigma = new double[n * n];

	// [U,Sigma,T]=svd(V,'econ');
	SvdUtility::getSVD(v, j, k, u, sigmas, tTransposed, sigma);
	SvdUtility::printMatrix("U", u, j, n);
	SvdUtility::printMatrix("Sigma", sigma, n, n);
	SvdUtility::printMatrix("T'", tTransposed, n, k);

	// NormS=norm(sigmas,2);
	double normS = SvdUtility::getEuclideanNorm(sigmas, n, n);
	std::cout << "NormS = " << normS << endl;

	// spatial complexity kk
	int kk = 0;
	for (int i = 0; i < n; ++i)
	{
		if (SvdUtility::getEuclideanNorm(sigmas, n, n - i) / normS > varepsilon1)
		{
			kk += 1;
		}
	}

	// Sigma(1:kk,1:kk)
	double* sigmaOneToN = new double[kk * kk];
	SvdUtility::resizeMatrix(sigma, sigmaOneToN, n, kk, 0, kk - 1);
	SvdUtility::printMatrix("Sigma(1:N,1:N)", sigmaOneToN, kk, kk);

	// T(:,1:n)'
	double* tOneToNTransposed = new double[kk * k];
	SvdUtility::resizeMatrix(tTransposed, tOneToNTransposed, n, kk, 0, k - 1);
	SvdUtility::printMatrix("T(:,1:N)^T", tOneToNTransposed, kk, k);

	// [N,K]=size(hatT);
	n = kk;
	std::cout << "N = rows = " << kk << endl << "K = cols = " << k << endl << endl;

	// hatT=Sigma(1:n,1:n)*T(:,1:n)'
	double* hatT = new double[n * k];
	SvdUtility::getMatrixProduct(sigmaOneToN, tOneToNTransposed, hatT, n, n, k, false);
	SvdUtility::printMatrix("hatT", hatT, n, k);

	int min = std::min(n, k - 1);

	// N > K => only the first N cols of hatU1 are computed
	double* hatU1 = new double[n * min];

	// N < K => only the first K rows of hatU2' are computed
	double* hatU2Transposed = new double[min * (k - 1)];

	double* hatSigmas = new double[min];
	double* hatSigma = new double[min * min];

	// hatT(:, 1 : K - 1)
	double* hatT0 = new double[n * (k - 1)];
	SvdUtility::resizeMatrix(hatT, hatT0, n, n, 0, k - 2);
	SvdUtility::printMatrix("hatT(:,1..K-1)", hatT0, n, k - 1);

	// [hatU1,hatSigma,hatU2]=svd(hatT(:,1:K-1),'econ');
	SvdUtility::getSVD(hatT0, n, k - 1, hatU1, hatSigmas, hatU2Transposed, hatSigma);
	SvdUtility::printMatrix("hatU1", hatU1, n, min);
	SvdUtility::printMatrix("hatSigma", hatSigma, min, min);
	SvdUtility::printMatrix("hatU2Transposed", hatU2Transposed, min, k - 1);

	// hatT(:,2:K)
	double* hatT1 = new double[n * (k - 1)];
	SvdUtility::resizeMatrix(hatT, hatT1, n, n, 1, k - 1);
	SvdUtility::printMatrix("hatT(:,2:K)", hatT1, n, k - 1);

	// hatU2	
	double* hatU2 = new double[(k - 1) * min];
	SvdUtility::transposeMatrix(hatU2Transposed, hatU2, min, k - 1);
	SvdUtility::printMatrix("hatU2", hatU2, k - 1, min);

	// hatT(:,2:K)*hatU2
	double* xHatU2 = new double[n * min];
	SvdUtility::getMatrixProduct(hatT1, hatU2, xHatU2, n, k - 1, min, false);
	SvdUtility::printMatrix("xHatU2", xHatU2, n, min);

	// inv(hatSigma)
	SvdUtility::getMatrixInverse(hatSigma, min);
	SvdUtility::printMatrix("hatSigma^-1", hatSigma, min, min);

	// hatT(:,2:K)*hatU2*inv(hatSigma)
	double* xHatSigmaInverse = new double[n * min];
	SvdUtility::getMatrixProduct(xHatU2, hatSigma, xHatSigmaInverse, n, min, min, false);
	SvdUtility::printMatrix("xHatSigmaInverse", xHatSigmaInverse, n, min);

	// hatU1'
	double* hatU1Transposed = new double[min * n];
	SvdUtility::transposeMatrix(hatU1, hatU1Transposed, n, min);
	SvdUtility::printMatrix("hatU1^T", hatU1Transposed, min, n);

	// hatR=hatT(:,2:K)*hatU2*inv(hatSigma)*hatU1'
	double* hatR = new double[n * n];
	SvdUtility::getMatrixProduct(xHatSigmaInverse, hatU1Transposed, hatR, n, min, n, false);
	SvdUtility::printMatrix("hatR", hatR, n, n);

	// double array hatR to complex double array hatRComplex
	double _Complex* hatRComplex = new double _Complex[n * n];
	for (int col = 0; col < n; ++col)
	{
		for (int row = 0; row < n; ++row)
		{
			hatRComplex[col * n + row] = hatR[col * n + row];
		}
	}

	// [Q,MM]=eig(hatR)
	// eigenvalues=diag(MM)
	double _Complex* eigenvalues = new double _Complex[n];
	double _Complex* q = new double _Complex[n * n];
	double _Complex* mm = new double _Complex[n * n];
	for (int i = 0; i < n; ++i)
	{
		mm[i * n + i] = eigenvalues[i];
	}
	SvdUtility::getEigen(hatRComplex, n, eigenvalues, q);
	SvdUtility::printMatrix("eigenvalues", eigenvalues, n, 1);
	SvdUtility::printMatrix("Q", q, n, n);
	//SvdUtility::printMatrix("MM", mm, n, n);

	// M=length(eigenvalues)
	//int m = n;

	// qq=log(eigenvalues)
	// deltas=real(qq)/Deltat
	// omegas=imag(qq)/Deltat
	double _Complex* qq = new double _Complex[n];
	double* deltas = new double[n];
	double* omegas = new double[n];
	for (int i = 0; i < n; ++i)
	{
		qq[i] = ::clog(eigenvalues[i]);
		deltas[i] = creal(qq[i]) / deltat;
		omegas[i] = cimag(qq[i]) / deltat;
	}
	SvdUtility::printMatrix("qq", qq, n, 1);
	SvdUtility::printMatrix("deltas", deltas, n, 1);
	SvdUtility::printMatrix("omegas", omegas, n, 1);

	/*double* mmm = new double[m * k * m];
	double* bb = new double[m * k];

	for (int i = 0; i < k, ++i)
	{

	}

	SvdUtility::printMatrix("Mm", mmm, m * k, m);
	SvdUtility::printMatrix("Bb", bb, m * k, 1);*/

	

			/*
			for(std::vector<double>::iterator it = result.begin(); it!=result.end(); ++it)
			{
			  std::cout << ' ' << *it;
			}

			std::cout << "/n" << std::endl;


			std::cout << "size Sigma: " << s.size() << std::endl;
			std::cout << "size U: " << u.size() << std::endl;
			std::cout << "size V transposed: " << vt.size() << std::endl;
			SvdUtility::writeCSV("./out/SVDresult.csv", result, columnCount, columnCount);
			*/

	return EXIT_SUCCESS;
}
