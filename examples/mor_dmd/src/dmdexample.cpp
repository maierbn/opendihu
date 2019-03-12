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
	double varepsilon = 0.0;
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
	int j = 6;
	int k = 5;
	double v[] {
			  8.79,  6.11, -9.15,  9.57, -3.49,  9.84,
			  9.93,  6.91, -7.93,  1.64,  4.02,  0.15,
			  9.83,  5.04,  4.86,  8.83,  9.80, -8.99,
			  5.45, -0.27,  4.85,  0.74, 10.00, -6.02,
			  3.16,  7.98,  3.01,  5.80,  4.27, -5.31
	};

	// complex square matrix
	/*int j = 5;
	int k = 5;
	double v[]{
			  8.79,  6.11, -9.15,  9.57, -3.49,
			  9.93,  6.91, -7.93,  1.64,  4.02,
			  9.83,  5.04,  4.86,  8.83,  9.80,
			  5.45, -0.27,  4.85,  0.74, 10.00,
			  3.16,  7.98,  3.01,  5.80,  4.27
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
	std::cout << "J = rows = " << j << endl << "K = cols =  " << k << endl << endl;
	// SvdUtility::printMatrix("V", v, j, k);

	// %%%%%%%%%%%%%%%%%%%%%%%%%
	// [U,Sigma,T]=svd(V,'econ')
	// %%%%%%%%%%%%%%%%%%%%%%%%%

	// n=length(sigmas);
	int n = std::min(j, k);

	// J > K => only the first K cols of U are computed
	double* u = new double[j * n];

	// J < K => only the first J rows of T^T are computed
	double* tTransposed = new double[n * k];

	double* sigmas = new double[n];
	double* sigma = new double[n * n];

	// [U,Sigma,T]=svd(V,'econ')
	SvdUtility::getSVD(v, j, k, u, sigmas, tTransposed, sigma);
	// SvdUtility::printMatrix("U", u, j, n);
	// SvdUtility::printMatrix("Sigma", sigma, n, n);
	// SvdUtility::printMatrix("T'", tTransposed, n, k);

	// %%%%%%%%%%%%%%%%%%%%%%
	// Spatial complexity: kk
	// %%%%%%%%%%%%%%%%%%%%%%

	// NormS=norm(sigmas,2);
	double normS = SvdUtility::getEuclideanNorm(sigmas, n, n);
	std::cout << "NormS = " << normS << endl;

	int kk = 0;
	for (int i = 0; i < n; ++i)
	{
		if (SvdUtility::getEuclideanNorm(sigmas, n, n - i) / normS > varepsilon1)
		{
			kk += 1;
		}
	}
	std::cout << "spatial complexity kk = " << kk << endl << endl;

	// U=U(:,1:kk)
	double* uOneTokk = new double[j * kk];
	SvdUtility::resizeMatrix(u, uOneTokk, j, j, 0, kk - 1);

	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// Create reduced snapshots matrix hatT
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	// Sigma(1:kk,1:kk)
	double* sigmaOneToN = new double[kk * kk];
	SvdUtility::resizeMatrix(sigma, sigmaOneToN, n, kk, 0, kk - 1);
	// SvdUtility::printMatrix("Sigma(1:kk,1:kk)", sigmaOneToN, kk, kk);

	// T(:,1:kk)'
	double* tOneToNTransposed = new double[kk * k];
	SvdUtility::resizeMatrix(tTransposed, tOneToNTransposed, n, kk, 0, k - 1);
	// SvdUtility::printMatrix("T(:,1:kk)^T", tOneToNTransposed, kk, k);

	// [N,K]=size(hatT);
	n = kk;
	std::cout << "N = rows = " << kk << endl << "K = cols = " << k << endl << endl;

	// hatT=Sigma(1:n,1:n)*T(:,1:n)'
	double* hatT = new double[n * k];
	SvdUtility::getMatrixProduct(sigmaOneToN, tOneToNTransposed, hatT, n, n, k, false);
	// SvdUtility::printMatrix("hatT", hatT, n, k);

	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// [hatU1,hatSigma,hatU2]=svd(hatT(:,1:K-1),'econ')
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
	// SvdUtility::printMatrix("hatT(:,1..K-1)", hatT0, n, k - 1);

	// [hatU1,hatSigma,hatU2]=svd(hatT(:,1:K-1),'econ')
	SvdUtility::getSVD(hatT0, n, k - 1, hatU1, hatSigmas, hatU2Transposed, hatSigma);
	// SvdUtility::printMatrix("hatU1", hatU1, n, min);
	// SvdUtility::printMatrix("hatSigma", hatSigma, min, min);
	// SvdUtility::printMatrix("hatU2Transposed", hatU2Transposed, min, k - 1);

	// %%%%%%%%%%%%%%%%%%%%%%%%%%
	// Calculate Koopman operator
	// %%%%%%%%%%%%%%%%%%%%%%%%%%

	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// hatR=hatT(:,2:K)*hatU2*inv(hatSigma)*hatU1'
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	// hatT(:,2:K)
	double* hatT1 = new double[n * (k - 1)];
	SvdUtility::resizeMatrix(hatT, hatT1, n, n, 1, k - 1);
	// SvdUtility::printMatrix("hatT(:,2:K)", hatT1, n, k - 1);

	// hatU2	
	double* hatU2 = new double[(k - 1) * min];
	SvdUtility::transposeMatrix(hatU2Transposed, hatU2, min, k - 1);
	// SvdUtility::printMatrix("hatU2", hatU2, k - 1, min);

	// hatT(:,2:K)*hatU2
	double* xHatU2 = new double[n * min];
	SvdUtility::getMatrixProduct(hatT1, hatU2, xHatU2, n, k - 1, min, false);
	// SvdUtility::printMatrix("xHatU2", xHatU2, n, min);

	// inv(hatSigma)
	SvdUtility::getMatrixInverse(hatSigma, min);
	// SvdUtility::printMatrix("hatSigma^-1", hatSigma, min, min);

	// hatT(:,2:K)*hatU2*inv(hatSigma)
	double* xHatSigmaInverse = new double[n * min];
	SvdUtility::getMatrixProduct(xHatU2, hatSigma, xHatSigmaInverse, n, min, min, false);
	// SvdUtility::printMatrix("xHatSigmaInverse", xHatSigmaInverse, n, min);

	// hatU1'
	double* hatU1Transposed = new double[min * n];
	SvdUtility::transposeMatrix(hatU1, hatU1Transposed, n, min);
	// SvdUtility::printMatrix("hatU1^T", hatU1Transposed, min, n);

	// hatR=hatT(:,2:K)*hatU2*inv(hatSigma)*hatU1'
	double* hatR = new double[n * n];
	SvdUtility::getMatrixProduct(xHatSigmaInverse, hatU1Transposed, hatR, n, min, n, false);
	// SvdUtility::printMatrix("hatR", hatR, n, n);

	// double array hatR to complex double array hatRComplex
	double _Complex* hatRComplex = new double _Complex[n * n];
	SvdUtility::doubleToComplex(hatR, hatRComplex, n, n);
	// SvdUtility::printMatrix("hatR", hatRComplex, n, n);

	// %%%%%%%%%%%%%%%%
	// [Q,MM]=eig(hatR)
	// %%%%%%%%%%%%%%%%

	double _Complex* eigenvalues = new double _Complex[n];
	double _Complex* q = new double _Complex[n * n];

	// [Q,MM]=eig(hatR)
	SvdUtility::getEigen(hatRComplex, n, eigenvalues, q);

	// eigenvalues=diag(MM)
	double _Complex* mm = new double _Complex[n * n];
	for (int i = 0; i < n; ++i)
	{
		mm[i * n + i] = eigenvalues[i];
	}
	// SvdUtility::printMatrix("eigenvalues", eigenvalues, n, 1);
	SvdUtility::printMatrix("Q", q, n, n);
	// SvdUtility::printMatrix("MM", mm, n, n);

	// M=length(eigenvalues)
	int m = n;

	// qq=log(eigenvalues)
	// deltas=real(qq)/Deltat
	// omegas=imag(qq)/Deltat
	double _Complex* qq = new double _Complex[n];
	double* deltas = new double[m];
	double* omegas = new double[m];
	for (int i = 0; i < m; ++i)
	{
		qq[i] = ::clog(eigenvalues[i]);
		deltas[i] = creal(qq[i]) / deltat;
		omegas[i] = cimag(qq[i]) / deltat;
	}
	// SvdUtility::printMatrix("qq", qq, n, 1);
	// SvdUtility::printMatrix("deltas", deltas, m, 1);
	// SvdUtility::printMatrix("omegas", omegas, m, 1);

	// %%%%%%%%%%%%%%%%%%%%
	// Calculate amplitudes
	// %%%%%%%%%%%%%%%%%%%%

	// Mm=zeros(M*K,M)
	double _Complex* mmm = new double _Complex[m * k * m];
	SvdUtility::setZero(mmm, m * k, m);

	// Bb=zeros(M*K,1)
	double _Complex* bb = new double _Complex[m * k];
	SvdUtility::setZero(bb, m * k, 1);

	// MM^(k-1)
	double _Complex* mmPower = new double _Complex[m * m];

	// Mm(1+(k-1)*M:k*M,:)
	double _Complex* mmmMxM = new double _Complex[m * m];

	for (int kay = 0; kay < k; ++kay)
	{
		SvdUtility::setZero(mmmMxM, m, m);
		SvdUtility::getMatrixPower(mm, mmPower, m, kay);
		// SvdUtility::printMatrix("MM^(k-1)", mmPower, m, m);
		// Q*(MM^(k-1))
		SvdUtility::getMatrixProduct(q, mmPower, mmmMxM, m, m, m);
		// SvdUtility::printMatrix("Mm(1+(k-1)*M:k*M,:)", mmmMxM, m, m);
		// Mm(1+(k-1)*M:k*M,:)=Q*(MM^(k-1))
		for (int row = 0; row < m; ++row)
		{
			for (int col = 0; col < m; ++col)
			{
				mmm[row + kay * m + col * k * m] = mmmMxM[row + col * m];
			}
			bb[row + kay * m] = hatT[row + kay * m];
		}
	}

	// SvdUtility::printMatrix("Mm", mmm, m * k, m);
	// SvdUtility::printMatrix("Bb", bb, m * k, 1);

	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// [Ur,Sigmar,Vr]=svd(Mm,'econ')
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	double _Complex* ur = new double _Complex[m * k * m];
	double* sigmars = new double[m];
	double _Complex* vrTransposed = new double _Complex[m * m];
	double* sigmar = new double[m * m];

	// [Ur,Sigmar,Vr]=svd(Mm,'econ')
	SvdUtility::getSVD(mmm, m * k, m, ur, sigmars, vrTransposed, sigmar);
	// SvdUtility::printMatrix("Ur", ur, m * k, m);
	// SvdUtility::printMatrix("Sigmar", sigmar, m, m);
	// SvdUtility::printMatrix("VrTransposed", vrTransposed, m, m);

	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// a=Vr*mldivide(Sigmar,Ur'*Bb)
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%

	// Ur'
	double _Complex* urTransposed = new double _Complex[m * m * k];
	SvdUtility::transposeMatrix(ur, urTransposed, m * k, m);

	// Ur'*Bb
	double _Complex* xBb = new double _Complex[m];
	SvdUtility::getMatrixProduct(urTransposed, bb, xBb, m, m * k, 1);
	// SvdUtility::printMatrix("Ur'*Bb", xBb, m, 1);

	// Sigmar\(Ur'*Bb)
	double _Complex* sigmarComplex = new double _Complex[m * m];
	SvdUtility::doubleToComplex(sigmar, sigmarComplex, m, m);
	SvdUtility::getMatrixLeftDivision(sigmarComplex, xBb, m, 1);
	// SvdUtility::printMatrix("mldivide(Sigmar,Ur'*Bb)", xBb, m, 1);

	//Vr
	double _Complex* vr = new double _Complex[m * m];
	SvdUtility::transposeMatrix(vrTransposed, vr, m, m);
	// SvdUtility::printMatrix("Vr", vr, m, m);

	// a=Vr*mldivide(Sigmar,Ur'*Bb)
	double _Complex* a = new double _Complex[m];
	SvdUtility::setZero(a, m, 1);
	SvdUtility::getMatrixProduct(vr, xBb, a, m, m, 1);
	SvdUtility::printMatrix("a", a, m, 1);

	// u=zeros(M:M)
	double _Complex* uComplex = new double _Complex[m * m];
	SvdUtility::setZero(uComplex, m, m);

	// u(:,m)=a(m)*Q(:,m)
	for (int col = 0; col < m; ++col)
	{
		for (int row = 0; row < m; ++row)
		{
			uComplex[row + col * m] = a[col] * q[row + col * m];
			// std::cout << creal(a[col] * q[row + col * m]) << " + " << cimag(a[col] * q[row + col * m]) << "   ";
		}
		// std::cout << endl;
	}
	SvdUtility::printMatrix("u", uComplex, m, m);

	// amplitude=zeros(M,1)
	double _Complex* amplitude = new double _Complex[m];
	SvdUtility::setZero(amplitude, m, 1);

	// aca=U*u(:,m)
	// amplitude(m)=norm(aca(:),2)/sqrt(J)
	double _Complex* uOneTokkComplex = new double _Complex[j * kk];
	// SvdUtility::printMatrix("U", uOneTokkComplex, j, kk);
	SvdUtility::doubleToComplex(uOneTokk, uOneTokkComplex, j, kk);
	double _Complex* aca = new double _Complex[j];
	double _Complex* uColm = new double _Complex[m];
	for (int col = 0; col < m; ++col)
	{
		SvdUtility::resizeMatrix(uComplex, uColm, m, m, col, col);
		// SvdUtility::printMatrix("u(:,m)", uColm, m, 1);
		SvdUtility::setZero(aca, j, 1);
		SvdUtility::getMatrixProduct(uOneTokkComplex, uColm, aca, j, m, 1);
		// SvdUtility::printMatrix("aca", aca, j, 1);
		amplitude[col] = SvdUtility::getEuclideanNorm(aca, j, j) / sqrt(j);
	}
	// SvdUtility::printMatrix("amplitude", amplitude, m, 1);

	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// UU=[u;deltas';omegas';amplitude']'
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	// transpose u
	double _Complex* uComplexTransposed = new double _Complex[m * m];
	SvdUtility::transposeMatrix(uComplex, uComplexTransposed, m, m);
	// SvdUtility::printMatrix("uTransposed", uComplexTransposed, m, m);

	double _Complex* deltasComplex = new double _Complex[m];
	SvdUtility::doubleToComplex(deltas, deltasComplex, m, 1);
	// SvdUtility::printMatrix("deltas", deltasComplex, m, 1);
	double _Complex* omegasComplex = new double _Complex[m];
	SvdUtility::doubleToComplex(omegas, omegasComplex, m, 1);
	// SvdUtility::printMatrix("omegas", omegasComplex, m, 1);

	// UU=[u;deltas';omegas';amplitude']'
	double _Complex* uu = new double _Complex[m * (m + 3)];
	SvdUtility::concatenateMatrices(uComplexTransposed, deltasComplex, uu, m, m);
	SvdUtility::concatenateMatrices(uu, omegasComplex, uu, m, m + 1);
	SvdUtility::concatenateMatrices(uu, amplitude, uu, m, m + 2);
	// SvdUtility::printMatrix("UU", uu, m, m + 3);

	// UU1=sortrows(UU,-(M+3))
	double _Complex* uu1 = new double _Complex[m * (m + 3)];
	SvdUtility::sortMatrix(uu, uu1, m, m + 3);
	// SvdUtility::printMatrix("UU1", uu1, m, m + 3);

	// UU=UU1'
	double _Complex* uu1Transposed = new double _Complex[(m + 3) * m];
	SvdUtility::transposeMatrix(uu1, uu1Transposed, m, m + 3);
	// SvdUtility::printMatrix("UU", uu1Transposed, m + 3, m);

	// u=UU(1:M,:)
	SvdUtility::resizeMatrix(uu1Transposed, uComplex, m + 3, m, 0, m - 1);
	// SvdUtility::printMatrix("u", uComplex, m, m);

	// deltas=UU(M+1,:)
	// omegas=UU(M+2,:)
	// amplitude=UU(M+3,:)
	for (int col = 0; col < m; ++col)
	{
		deltasComplex[col] = uu1Transposed[m + col * (m + 3)];
		omegasComplex[col] = uu1Transposed[m + 1 + col * (m + 3)];
		amplitude[col] = uu1Transposed[m + 2 + col * (m + 3)];
	}
	// SvdUtility::printMatrix("deltas", deltasComplex, m, 1);
	// SvdUtility::printMatrix("omegas", omegasComplex, m, 1);
	// SvdUtility::printMatrix("amplitude", amplitude, m, 1);

	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// Spectral complexity: number of DMD modes
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	int kk2 = 0;
	for (int em = 0; em < m; ++em)
	{
		if (creal(amplitude[em] / amplitude[0]) > varepsilon)
		{
			kk2 += 1;
		}
	}
	std::cout << "Spectral complexity: number of DMD modes kk2 = " << kk2 << endl << endl;

	// u=u(:,1:kk2)
	double _Complex* uOneTokk2 = new double _Complex[m * kk2];
	SvdUtility::resizeMatrix(uComplex, uOneTokk2, m, m, 0, kk2 - 1);
	// SvdUtility::printMatrix("u", uOneTokk2, m, kk2);

	// deltas=deltas(1:kk2)
	double _Complex* deltasOneTokk2 = new double _Complex[kk2];
	SvdUtility::resizeMatrix(deltasComplex, deltasOneTokk2, m, kk2, 0, 0);
	// SvdUtility::printMatrix("deltas", deltasOneTokk2, m, kk2);

	// omegas=omegas(1:kk2)
	double _Complex* omegasOneTokk2 = new double _Complex[kk2];
	SvdUtility::resizeMatrix(omegasComplex, omegasOneTokk2, m, kk2, 0, 0);
	// SvdUtility::printMatrix("omegas", omegasOneTokk2, m, kk2);

	// amplitude=amplitude(1:kk2)
	double _Complex* amplitudeOneTokk2 = new double _Complex[kk2];
	SvdUtility::resizeMatrix(amplitude, amplitudeOneTokk2, m, kk2, 0, 0);
	// SvdUtility::printMatrix("amplitude", amplitudeOneTokk2, m, kk2);

	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// DeltasOmegAmpl=[deltas',omegas',amplitude']
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	double _Complex* deltasOmegAmpl = new double _Complex[kk2 * 3];
	for (int row = 0; row < kk2; ++row)
	{
		deltasOmegAmpl[row] = deltasOneTokk2[row];
		deltasOmegAmpl[row + kk2] = omegasOneTokk2[row];
		deltasOmegAmpl[row + 2 * kk2] = amplitudeOneTokk2[row];
	}
	// SvdUtility::printMatrix("DeltasOmegAmpl", deltasOmegAmpl, kk2, 3);

	//hatTreconst=zeros(N,K)
	double _Complex* hatTreconst = new double _Complex[n * k];
	SvdUtility::setZero(hatTreconst, n, k);

	double _Complex* hatTreconstCol = new double _Complex[n];
	for (int col = 0; col < k; ++col)
	{
		// hatTreconst(:,k)=ContReconst_SIADS(Time(k),Time(1),u,deltas,omegas);
		SvdUtility::setZero(hatTreconstCol, n, 1);
		SvdUtility::contReconst(col * deltat, 0, uOneTokk2, deltasOneTokk2, omegasOneTokk2, m, kk2, hatTreconstCol);
		for (int row = 0; row < n; ++row)
		{
			hatTreconst[row + col * n] = hatTreconstCol[row];
		}
	}
	// SvdUtility::printMatrix("hatTreconst", hatTreconst, n, k);

	// Vreconst=U*hatTreconst
	double _Complex* vreconst = new double _Complex[j * k];
	SvdUtility::getMatrixProduct(uOneTokkComplex, hatTreconst, vreconst, j, n, k);
	SvdUtility::printMatrix("Vreconst", vreconst, j, k);

	// NormV=norm(V(:),2)
	double normV = sigmas[0];
	std::cout << "NormV = " << normV << endl;

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
