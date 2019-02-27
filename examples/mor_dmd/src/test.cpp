#include <iostream>
#include <cstdlib>
#include <vector>
#include "opendihu.h"
#include "utility/svd_utility.h"

int main(int argc, char *argv[])
{
	// input data
	//std::vector<double> parsedCsv = SvdUtility::readCSV("./out/data.csv");
	//double v[parsedCsv.size()];
	//copy(parsedCsv.begin(), parsedCsv.end(), v);

	//int rowCount = SvdUtility::getCSVRowCount("./out/data.csv");
	//int columnCount = SvdUtility::getCSVColumnCount("./out/data.csv");

	double varepsilon1 = 0;

	int j = 3;
	int k = 3;
	double v[j * k] = {
		1, -1,  0,
		0, -2,  1, 
		1,  0, -1
	};

	// [J,K]=size(V);
	/*int j = 6;
	int k = 5;
	// 5x6
	double v[j * k] = {
			  8.79,  6.11, -9.15,  9.57, -3.49,  9.84,
			  9.93,  6.91, -7.93,  1.64,  4.02,  0.15,
			  9.83,  5.04,  4.86,  8.83,  9.80, -8.99,
			  5.45, -0.27,  4.85,  0.74, 10.00, -6.02,
			  3.16,  7.98,  3.01,  5.80,  4.27, -5.31
	};*/

	SvdUtility::printMatrix("V", v, j, k);

	// 5x5
	/*int j = 5;
	int k = 5;
	double v[j * k] = {
			  8.79,  6.11, -9.15,  9.57, -3.49,
			  9.93,  6.91, -7.93,  1.64,  4.02,
			  9.83,  5.04,  4.86,  8.83,  9.80,
			  5.45, -0.27,  4.85,  0.74, 10.00,
			  3.16,  7.98,  3.01,  5.80,  4.27
	};*/

	// 6x5
	/*int j = 5;
	int k = 6;
	double v[j * k] = {
			  8.79,  6.11, -9.15,  9.57, -3.49,
			  9.93,  6.91, -7.93,  1.64,  4.02,
			  9.83,  5.04,  4.86,  8.83,  9.80,
			  5.45, -0.27,  4.85,  0.74, 10.00,
			  3.16,  7.98,  3.01,  5.80,  4.27,
			  9.84,  0.15, -8.99, -6.02, -5.31
	};*/

	// print input matrix sizes, j cols, k rows
	std::cout << "J = " << j << endl << "K = " << k << endl << endl;

	// n=length(sigmas);
	int n = std::min(j, k);

	// m > n � Only the first n cols of U are computed, and S is n - by - n.
	// m = n � svd(A, 'econ') is equivalent to svd(A).
	// m < n � Only the first m cols of V are computed, and S is m - by - m.
	int colsU = j;
	int colsTTransposed = k;

	if (j > k) {
		colsU = k;
	}
	else if (j < k) {
		colsTTransposed = j;
	}

	// [U,Sigma,T]=svd(V,'econ');
	double u[j * colsU];
	double sigmas[n];
	double tTransposed[k * colsTTransposed];
	double sigma[n * n];
	SvdUtility::getSVD(v, j, k, u, sigmas, tTransposed, sigma);

	// sigmas=diag(Sigma);
	/*double sigma[n * n] = {};
	for (int row = 0; row < n; ++row)
	{
		for (int column = 0; column < n; ++column)
		{
			if (row == column)
			{
				sigma[row + column * n] = sigmas[row];
			} 
		}
	}*/

	// print resulting matrices
	SvdUtility::printMatrix("U", u, j, colsU);
	SvdUtility::printMatrix("Sigma", sigma, n, n);
	SvdUtility::printMatrix("T^T", tTransposed, k, colsTTransposed);

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
	std::cout << "kk = " << kk << endl << endl;

	// hatT=Sigma(1:kk,1:kk)*T(:,1:kk)'
	double hatT[kk * colsTTransposed] = {};
	SvdUtility::resizeMatrix(sigma, sigma, n, n, 0, kk - 1, 0, kk - 1);
	SvdUtility::resizeMatrix(tTransposed, tTransposed, k, colsTTransposed, 0, kk - 1, 0, colsTTransposed - 1);
	SvdUtility::printMatrix("Sigma(1:kk,1:kk)", sigma, kk, kk);
	SvdUtility::printMatrix("T(..1:kk)", tTransposed, kk, colsTTransposed);
	SvdUtility::getMatrixProduct(sigma, tTransposed, hatT, kk, kk, colsTTransposed);
	SvdUtility::printMatrix("hatT", hatT, kk, colsTTransposed);

	// [N,K]=size(hatT);
	n = kk;
	k = colsTTransposed;
 
	// [hatU1,hatSigma,hatU2]=svd(hatT(:,1:K-1),'econ');
	int min = std::min(k - 1, n);
	int colsHatU1 = n;
	int colsHatU2Transposed = n;s
	/*if (k - 1 < n) {
		colsHatU1 = k - 1;
	}
	else if (k - 1 > n) {
		colsHatU2Transposed = n;
	}*/
	double hatU1[n * colsHatU1];
	double hatSigmas[min];
	double hatSigma[min * min];
	double hatU2Transposed[k - 1 * colsHatU2Transposed];
	double hatT0[n * k - 1];
	std::cout << "N = " << n << endl << "K = " << k << endl << endl;
	SvdUtility::resizeMatrix(hatT, hatT0, n, k, 0, n - 1, 0, k - 2);
	SvdUtility::printMatrix("hatT(:,1..K-1)", hatT0, n, k - 1);
	SvdUtility::getSVD(hatT0, k - 1, n, hatU1, hatSigmas, hatU2Transposed, hatSigma);
	SvdUtility::printMatrix("hatU1", hatU1, n, colsHatU1);
	SvdUtility::printMatrix("hatSigma", hatSigma, min, min);
	SvdUtility::printMatrix("hatU2Transposed", hatU2Transposed, k - 1, colsHatU2Transposed);

	// hatR=hatT(:,2:K)*hatU2*inv(hatSigma)*hatU1';
	double hatT1[n * k - 1];
	//double hatR[n * k];
	double hatU2[colsHatU2Transposed * k - 1];
	double hatT1xHatU2[n * n];
	SvdUtility::resizeMatrix(hatT, hatT1, n, k, 0, n - 1, 1, k - 1);
	SvdUtility::printMatrix("hatT(:,2:K)", hatT1, n, k - 1);
	SvdUtility::transposeMatrix(hatU2Transposed, hatU2, k - 1, colsHatU2Transposed);
	SvdUtility::printMatrix("hatU2", hatU2, colsHatU2Transposed, k - 1);
	SvdUtility::getMatrixProduct(hatT1, hatU2, hatT1xHatU2, n, colsHatU2Transposed, n);
	SvdUtility::printMatrix("hatT1xHatU2", hatT1xHatU2, n, n);
	SvdUtility::getMatrixInverse(hatSigma, min);
	SvdUtility::printMatrix("hatSigmaInverse", hatSigma, min, min);
	//SvdUtility::getMatrixProduct(hatT1xHatU2, hatSigmaInverted, xHatSigmaInverted, n, n, n);

	
	
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
