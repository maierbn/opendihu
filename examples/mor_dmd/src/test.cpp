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

	double varepsilon1 = 0.3;

	// [J,K]=size(V);
	int j = 6;
	int k = 5;
	// 5x6
	double v[j * k] = {
			  8.79,  6.11, -9.15,  9.57, -3.49,  9.84,
			  9.93,  6.91, -7.93,  1.64,  4.02,  0.15,
			  9.83,  5.04,  4.86,  8.83,  9.80, -8.99,
			  5.45, -0.27,  4.85,  0.74, 10.00, -6.02,
			  3.16,  7.98,  3.01,  5.80,  4.27, -5.31
	};

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

	// print input matrix sizes, j columns, k rows
	std::cout << "J = " << j << endl << "K = " << k << endl << endl;

	// n=length(sigmas);
	int n = std::min(j, k);

	// m > n — Only the first n columns of U are computed, and S is n - by - n.
	// m = n — svd(A, 'econ') is equivalent to svd(A).
	// m < n — Only the first m columns of V are computed, and S is m - by - m.
	int columnsU = j;
	int columnsTT = k;

	if (j > k) {
		columnsU = k;
	}
	else if (j < k) {
		columnsTT = j;
	}

	// [U,Sigma,T]=svd(V,'econ');
	double u[j * columnsU];
	double sigmas[n];
	double tT[k * columnsTT];
	SvdUtility::getSVD(v, j, k, u, sigmas, tT);

	// sigmas=diag(Sigma);
	double sigma[n * n] = {};
	for (int row = 0; row < n; ++row)
	{
		for (int column = 0; column < n; ++column)
		{
			if (row == column)
			{
				sigma[row + column * n] = sigmas[row];
			} 
		}
	}

	// print resulting matrices
	SvdUtility::printMatrix("U", u, j, columnsU);
	SvdUtility::printMatrix("Sigma", sigma, n, n);
	SvdUtility::printMatrix("T^T", tT, k, columnsTT);

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

	// hatT=Sigma(1:kk,1:kk)*T(..1:kk)'
	double hatT[kk * columnsTT] = {};
	SvdUtility::resizeMatrix(sigma, sigma, n, n, 0, kk - 1, 0, kk - 1);
	SvdUtility::resizeMatrix(tT, tT, k, columnsTT, 0, kk - 1, 0, columnsTT - 1);
	SvdUtility::printMatrix("Sigma(1:kk,1:kk)", sigma, kk, kk);
	SvdUtility::printMatrix("T(..1:kk)", tT, kk, columnsTT);
	SvdUtility::getMatrixProduct(sigma, tT, hatT, kk, kk, columnsTT);
	SvdUtility::printMatrix("hatT", hatT, kk, columnsTT);

	n = kk;
	k = columnsTT;
 
	int minKN = std::min(k, n);

	int columnsHatU1 = k;
	int columnsHatU2 = n;

	if (k > n) {
		columnsHatU1 = n;
	}
	else if (k < n) {
		columnsHatU2 = k;
	}
	double hatU1[k * columnsHatU1];
	double hatSigma[minKN];
	double hatU2[n * columnsHatU2];
	std::cout << "N = " << n << endl << "K = " << k << endl;
	SvdUtility::resizeMatrix(hatT, hatT, n, k, 0, n - 1, 0, k - 2);
	SvdUtility::getSVD(hatT, k, n, hatU1, hatSigma, hatU2);
	SvdUtility::printMatrix("hatSigma", hatSigma, minKN, 1);
	/**
					// N = kk, K = rowCount

					std::cout << "hatT0" << endl;
					double hatT0[kk * (rowCount - 1)];
					for (int i = 0; i < kk; ++i)
					{
						for (int j = 0; j < rowCount - 1; ++j)
						{
							hatT0[j + i * (rowCount - 1)] = hatT[j + i * rowCount];
							std::cout << hatT0[j + i * (rowCount - 1)] << " ";
						}
						std::cout << endl;
					}
					std::cout << endl;

					n = std::min(kk, rowCount - 1);
					double hatU1[(rowCount - 1) * (rowCount - 1)];
					double hatSigma[n];
					double hatU2Transposed[kk * kk];
					SvdUtility::getSVD(hatT0, rowCount - 1, kk, hatU1, hatSigma, hatU2Transposed);

					std::cout << "hatU1" << endl;
					for (int i = 0; i < rowCount - 1; ++i)
					{
						for (int j = 0; j < rowCount - 1; ++j)
						{
							std::cout << hatU1[j + i * rowCount] << " ";
						}
						std::cout << endl;
					}

					std::cout << endl;

					std::cout << "hatSigma" << endl;
					for (int i = 0; i < n; ++i)
					{
						std::cout << hatSigma[i] << " ";
					}
					std::cout << endl << endl;

					std::cout << "hatU2Transposed" << endl;
					for (int i = 0; i < kk; ++i)
					{
						for (int j = 0; j < kk; ++j)
						{
							std::cout << hatU2Transposed[j + i * kk] << " ";
						}
						std::cout << endl;
					}
			**/
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
