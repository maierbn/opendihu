#include "opendihu.h"
#include "utility/svd_utility.h"
#include "utility/dmd_utility.h"

#include <iostream>
#include <cstdlib>
#include <vector>
#include <complex>
#include <math.h>
//#include <chrono>

// using namespace std;

int main(int argc, char *argv[])
{
	// input data
	std::string inputData = "./out/snapshots.csv";

	std::vector<double> parsedCsv = SvdUtility::readCSV(inputData);
	int j = SvdUtility::getCSVColumnCount(inputData);
	int k = SvdUtility::getCSVRowCount(inputData);
	std::cout << "J = " << j << endl << "K = " << k << endl;
	double* v = new double[j * k];
	std::copy(parsedCsv.begin(), parsedCsv.end(), v);

	double deltat = 0.1;
	double epsilon1 = exp(-2);
	double epsilon0 = exp(-1);

	//DmdUtility::printMatrix("V", v, j, k);

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

	//DmdUtility::printMatrix("\\hat{T}", hatT, kSpat, k);

	// %%%%%%%%%%%%%%%%%%%%%%%%%%
	// Calculate Koopman operator
	// %%%%%%%%%%%%%%%%%%%%%%%%%%

	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// hatR=hatT(:,2:K)*hatU2*inv(hatSigma)*hatU1'
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	double* hatR = new double[kSpat * kSpat];
	DmdUtility::getKoopmanOperator(hatT, hatR, kSpat, k);
	//DmdUtility::printMatrix("\\hat{R}", hatR, kSpat, kSpat);

	// %%%%%%%%%%%%%%%%
	// [Q,MM]=eig(hatR)
	// %%%%%%%%%%%%%%%%

	double _Complex* eigenvalues = new double _Complex[kSpat];
	double _Complex* q = new double _Complex[kSpat * kSpat];

	// [Q,MM]=eig(hatR)
	DmdUtility::getEigen(hatR, kSpat, eigenvalues, q);

	//DmdUtility::printMatrix("\\Lambda", eigenvalues, kSpat, 0);
	//DmdUtility::printMatrix("Q", q, kSpat, kSpat);

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

	//DmdUtility::printMatrix("\\delta", deltas, kSpat, 1);
	//DmdUtility::printMatrix("\\omega", omegas, kSpat, 1);

	// %%%%%%%%%%%%%%%%%%%%
	// Calculate amplitudes
	// %%%%%%%%%%%%%%%%%%%%

	double _Complex* a = new double _Complex[kSpat];

	DmdUtility::getAmplitudes(hatT, eigenvalues, q, kSpat, k, a);

	// DmdUtility::printMatrix("a", a, kSpat, 1);

	// u=zeros(M:M)
	double _Complex* uComplex = new double _Complex[kSpat * kSpat];
	double* amplitude = new double[kSpat];

	// U=U(:,1:kk)
	double* uOneTokSpat = new double[j * kSpat];

	DmdUtility::getDmdModes(uComplex, a, q, j, kSpat, u, uOneTokSpat, amplitude);

	// DmdUtility::printMatrix("u", u, kSpat, kSpat);

	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// UU=[u;deltas';omegas';amplitude']'
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// Spectral complexity: number of DMD modes
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	int kSpec = DmdUtility::getSpecComp(uComplex, deltas, omegas, amplitude, kSpat, epsilon0);

	// DmdUtility::printMatrix("\\alpha", amplitude, kSpat, 1);

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
	double normV = DmdUtility::getNorm(v, j, k);

	double* diff = new double[j * k];
	DmdUtility::setZero(diff, j, k);
	for (int row = 0; row < j; ++row)
	{
		for (int col = 0; col < k; ++col)
		{
			diff[row + col * j] = v[row + col * j] - vreconstReal[row + col * j];
		}
	}

	double rms = DmdUtility::getNorm(diff, j, k) / normV;
	std::cout << "RelativeerrorRMS = " << rms << endl;

	double max = DmdUtility::getInfNorm(diff, j, k) / DmdUtility::getInfNorm(v, j, k);
	std::cout << "RelativeerrorMax = " << max << endl;

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
