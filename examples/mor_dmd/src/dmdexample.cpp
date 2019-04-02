#include "opendihu.h"
#include "utility/svd_utility.h"
#include "utility/dmd_utility.h"

#include <iostream>
#include <cstdlib>
#include <vector>
#include <complex>
#include <math.h>

// using namespace std;

int main(int argc, char *argv[])
{
	// input data
	std::vector<double> parsedCsv = SvdUtility::readCSV("./out/data.csv");
	//int rowCount = DmdUtility::getCSVRowCount("data.csv");
	//std::cout << "rowCount = " << rowCount << endl;
	int j = SvdUtility::getCSVColumnCount("./out/data.csv");
	int k = SvdUtility::getCSVRowCount("./out/data.csv");
	std::cout << "J = " << j << endl << "K = " << k << endl;
	double* v = new double[j * k];
	std::copy(parsedCsv.begin(), parsedCsv.end(), v);
	//DmdUtility::readCSV(".out/data.csv", v, j, k);

	double epsilon1 = 0.1;
	double epsilon0 = 0.1;
	double deltat = 0.1;

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

	// complex wide and short matrix
	/*int j = 5;
	int k = 6;
	double v[]{
			  8.79,  6.11, -9.15,  9.57, -3.49,
			  9.93,  6.91, -7.93,  1.64,  4.02,
			  9.83,  5.04,  4.86,  8.83,  9.80,
			  5.45, -0.27,  4.85,  0.74, 10.00,
			  3.16,  7.98,  3.01,  5.80,  4.27,
			  9.84,  0.15, -8.99, -6.02, -5.31
	};*/

	// print input matrix sizes, j rows, k cols
	// std::cout << "J = rows = " << j << endl << "K = cols =  " << k << endl << endl;
	DmdUtility::printMatrix("V", v, j, k);

	// %%%%%%%%%%%%%%%%%%%%%%%%%
	// [U,Sigma,T]=svd(V,'econ')
	// %%%%%%%%%%%%%%%%%%%%%%%%%

	// n=length(sigmas);
	int n = std::min(j, k);

	// J > K => only the first K cols of U are computed
	double* u = new double[j * n];

	// J < K => only the first J rows of T^T are computed
	double* tT = new double[n * k];

	//double* sigmas = new double[n];
	double* sigma = new double[n * n];

	int kSpat = DmdUtility::getSpatComp(v, j, k, u, sigma, tT, epsilon1);

	// %%%%%%%%%%%%%%%%%%%%%%
	// Spatial complexity: kk
	// %%%%%%%%%%%%%%%%%%%%%%

	std::cout << "Spatial complexity: kSpat = " << kSpat << endl << endl;

	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// Create reduced snapshots matrix hatT
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	double* hatT = new double[kSpat * k];
	DmdUtility::getReducedSnapshotsMatrix(sigma, tT, hatT, n, k, kSpat);

	// %%%%%%%%%%%%%%%%%%%%%%%%%%
	// Calculate Koopman operator
	// %%%%%%%%%%%%%%%%%%%%%%%%%%

	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// hatR=hatT(:,2:K)*hatU2*inv(hatSigma)*hatU1'
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	double* hatR = new double[kSpat * kSpat];
	DmdUtility::getKoopmanOperator(hatT, hatR, kSpat, k);

	// %%%%%%%%%%%%%%%%
	// [Q,MM]=eig(hatR)
	// %%%%%%%%%%%%%%%%

	double _Complex* eigenvalues = new double _Complex[kSpat];
	double _Complex* q = new double _Complex[kSpat * kSpat];

	// [Q,MM]=eig(hatR)
	DmdUtility::getEigen(hatR, kSpat, eigenvalues, q);

	// eigenvalues=diag(MM)
	double _Complex* mm = new double _Complex[kSpat * kSpat];
	for (int i = 0; i < kSpat; ++i)
	{
		mm[i * kSpat + i] = eigenvalues[i];
	}

	// M=length(eigenvalues)
	//int kSpat = n;

	double* deltas = new double[kSpat];
	double* omegas = new double[kSpat];

	DmdUtility::getDeltaOmega(eigenvalues, deltas, omegas, kSpat, deltat);

	// %%%%%%%%%%%%%%%%%%%%
	// Calculate amplitudes
	// %%%%%%%%%%%%%%%%%%%%

	double _Complex* a = new double _Complex[kSpat];

	DmdUtility::getAmplitudes(hatT, eigenvalues, q, kSpat, k, a);

	// u=zeros(M:M)
	double _Complex* uComplex = new double _Complex[kSpat * kSpat];
	double* amplitude = new double[kSpat];

	// U=U(:,1:kk)
	double* uOneTokSpat = new double[j * kSpat];

	DmdUtility::getDmdModes(uComplex, a, q, j, kSpat, u, uOneTokSpat, amplitude);

	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// UU=[u;deltas';omegas';amplitude']'
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// Spectral complexity: number of DMD modes
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	int kSpec = DmdUtility::getSpecComp(uComplex, deltas, omegas, amplitude, kSpat, epsilon0);

	std::cout << "Spectral complexity: number of DMD modes kSpec = " << kSpec << endl << endl;
	
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// DeltasOmegAmpl=[deltas',omegas',amplitude']
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	double _Complex* uReduced = new double _Complex[kSpat * kSpec];
	double* deltasReduced = new double[kSpec];
	double* omegasReduced = new double[kSpec];

	DmdUtility::getDmdModesGrowthRatesFrequencies(uComplex, deltas, omegas, uReduced, deltasReduced, omegasReduced, kSpat, kSpat, kSpec);

	double* vreconstReal = new double[j * k];
	
	DmdUtility::reconstructSnapshots(uReduced, deltasReduced, omegasReduced, uOneTokSpat, vreconstReal, j, k, kSpec, kSpat, deltat);

	DmdUtility::printMatrix("Vreconst", vreconstReal, j, k);

	// NormV=norm(V(:),2)
	//double _Complex normV = SvdUtility::getMatrix2Norm(vreconstReal, j, k);
	//std::cout << "NormV = " << creal(normV) << endl;

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
