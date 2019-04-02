#include "utility/dmd_utility.h"
#include "utility/svd_utility.h"
#include <cblas.h>

#include <iostream>
#include <string>
#include <fstream>

int DmdUtility::getSpatComp(double input[], int rows, int cols, double leftSingVec[], double singVal[], double rightSingVec[], double epsilon)
{
	SvdUtility::getSVD(input, rows, cols, leftSingVec, singVal, rightSingVec);
	int min = std::min(rows, cols);
	double normS = DmdUtility::getL2norm(singVal, min, min);
	int spatComp = 0;
	for (int i = 0; i < min; ++i)
	{
		if (DmdUtility::getL2norm(singVal, min, min - i) / normS > epsilon)
		{
			spatComp += 1;
		}
	}
	return spatComp;
}

// takes real vector input with size entries as double[]
// returns the L2 norm of vector input over the specified range of entries as double
// sqrt(sum from k = size - range to size - 1 of v[k]^2)
double DmdUtility::getL2norm(double input[], int order, int range)
{
	double sum = 0;
	for (int k = order - range; k < order; ++k)
	{
		sum += input[order * k + k] * input[order * k + k];
	}
	return sqrt(sum);
}

// takes complex vector input with size entries as double _Complex[]
// returns the L2 norm of vector input over the specified range of entries as double
// sqrt(sum from k = size - range to size - 1 of Re(v[k])^2 + Im(v[k])^2)
double DmdUtility::getL2norm(double _Complex input[], int order, int range)
{
	double _Complex sum = 0;
	for (int k = order - range; k < order; ++k)
	{
		sum += creal(input[k]) * creal(input[k]) + cimag(input[k]) * cimag(input[k]);
	}
	return creal(csqrt(sum));
}

void DmdUtility::getReducedSnapshotsMatrix(double sigma[], double rightSingVec[], double hatT[], int min, int cols, int spatComp)
{
	double* sigmaResized = new double[spatComp * spatComp];
	DmdUtility::resizeMatrix(sigma, sigmaResized, min, spatComp, 0, spatComp - 1);

	double* rightSingVecResized = new double[spatComp * cols];
	DmdUtility::resizeMatrix(rightSingVec, rightSingVecResized, min, spatComp, 0, cols - 1);

	DmdUtility::getMatrixMult(sigmaResized, rightSingVecResized, hatT, spatComp, spatComp, cols);
}

// takes real matrix input with oldRows rows as double[]
// stores the entries of input in real matrix output with newRows rows and lastCol - firstCol + 1 columns as double[] starting with input[0, firstCol] and ending with input[newRows - 1, lastCol]
void DmdUtility::resizeMatrix(double input[], double output[], int oldRows, int newRows, int firstCol, int lastCol)
{
	for (int col = 0; col < lastCol - firstCol + 1; ++col)
	{
		for (int row = 0; row < newRows; ++row)
		{
			output[col * newRows + row] = input[(firstCol + col) * oldRows + row];
		}
	}
}

// takes complex matrix input with oldRows rows as double _Complex[]
// stores the entries of input in complex matrix output with newRows rows and lastCol - firstCol + 1 columns as double _Complex[] starting with input[0, firstCol] and ending with input[newRows - 1, lastCol]
void DmdUtility::resizeMatrix(double _Complex input[], double _Complex output[], int oldRows, int newRows, int firstCol, int lastCol)
{
	for (int col = 0; col < lastCol - firstCol + 1; ++col)
	{
		for (int row = 0; row < newRows; ++row)
		{
			output[col * newRows + row] = input[(firstCol + col) * oldRows + row];
		}
	}
}

// takes two real matrices inputA (rowsA x colsA_rowsB) and inputB (colsA_rowsB x colsB) each as double[]
// performs matrix multiplication utilizing cblas_dgemm
// stores if (trans) {inputA * inputB^T} else {inputA * inputB} in output (rowsA x colsB) as double[]
void DmdUtility::getMatrixMult(double inputA[], double inputB[], double output[], int rowsA, int colsA_rowsB, int colsB)
{
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, rowsA, colsB, colsA_rowsB, 1, inputA, rowsA, inputB, colsA_rowsB, 0, output, rowsA);
}

// takes two complex matrices inputA (rowsA x colsA_rowsB) and inputB (colsA_rowsB x colsB) each as double _Complex[]
// stores inputA * inputB in output (rowsA x colsB) as double _Complex[]
void DmdUtility::getMatrixMult(double _Complex inputA[], double _Complex inputB[], double _Complex output[], int rowsA, int colsA_rowsB, int colsB)
{
	// cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, rowsA, colsB, colsA_rowsB, 1, inputA, rowsA, inputB, colsA_rowsB, 0, output, rowsA);

	for (int col = 0; col < colsB; ++col)
	{
		for (int row = 0; row < rowsA; ++row)
		{
			output[row + col * rowsA] = 0;
			for (int cell = 0; cell < colsA_rowsB; ++cell)
			{
				output[row + col * rowsA] += inputA[row + cell * rowsA] * inputB[col * colsA_rowsB + cell];
			}
		}
	}
}

// takes one real matrix inputA (rowsA x colsA_rowsB) as double[] and one complex matrix inputB (colsA_rowsB x colsB) as double _Complex[]
// stores inputA * inputB in output (rowsA x colsB) as double _Complex[]
void DmdUtility::getMatrixMult(double inputA[], double _Complex inputB[], double _Complex output[], int rowsA, int colsA_rowsB, int colsB)
{
	// cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, rowsA, colsB, colsA_rowsB, 1, inputA, rowsA, inputB, colsA_rowsB, 0, output, rowsA);

	for (int col = 0; col < colsB; ++col)
	{
		for (int row = 0; row < rowsA; ++row)
		{
			output[row + col * rowsA] = 0;
			for (int cell = 0; cell < colsA_rowsB; ++cell)
			{
				output[row + col * rowsA] += inputA[row + cell * rowsA] * inputB[col * colsA_rowsB + cell];
			}
		}
	}
}

// takes two complex matrices inputA (rowsA x colsA_rowsB) and inputB (colsA_rowsB x colsB) each as double _Complex[]
// stores inputA * inputB in output (rowsA x colsB) as double[]
void DmdUtility::getMatrixMult(double _Complex inputA[], double _Complex inputB[], double output[], int rowsA, int colsA_rowsB, int colsB)
{
	// cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, rowsA, colsB, colsA_rowsB, 1, inputA, rowsA, inputB, colsA_rowsB, 0, output, rowsA);

	for (int col = 0; col < colsB; ++col)
	{
		for (int row = 0; row < rowsA; ++row)
		{
			output[row + col * rowsA] = 0;
			for (int cell = 0; cell < colsA_rowsB; ++cell)
			{
				output[row + col * rowsA] += creal(inputA[row + cell * rowsA] * inputB[col * colsA_rowsB + cell]);
			}
		}
	}
}

// takes real matrix input (rows x cols) as double[]
// prints name and input entry-wise
void DmdUtility::printMatrix(std::string name, double input[], int rows, int cols)
{
	cout << name << endl;
	for (int row = 0; row < rows; ++row)
	{
		for (int col = 0; col < cols; ++col)
		{
			if (input[col * rows + row] >= 0)
			{
				cout << " ";
			}
			cout << input[col * rows + row];
			if (col < cols - 1)
			{
				cout << "   ";
			}
		}
		cout << endl;
	}
	cout << endl;
}

// takes complex matrix input (rows x cols) as double _Complex[]
// prints name and input entry-wise
void DmdUtility::printMatrix(std::string name, double _Complex input[], int rows, int cols)
{
	cout << name << endl;
	for (int row = 0; row < rows; ++row)
	{
		for (int col = 0; col < cols; ++col)
		{
			if (creal(input[col * rows + row]) >= 0)
			{
				cout << " ";
			}
			cout << creal(input[col * rows + row]);
			if (cimag(input[col * rows + row]) < 0)
			{
				cout << " - " << -cimag(input[col * rows + row]) << "i";
			}
			else
			{
				cout << " + " << cimag(input[col * rows + row]) << "i";
			}
			if (col < cols - 1)
			{
				cout << "   ";
			}
		}
		cout << endl;
	}
	cout << endl;
}

// takes real matrix input (rows x cols) as double[]
// computes pseudoninverse of input by SVD
// computes real matrix output (rows x rows) where input(:,2:cols) = output * input(:,1:cols-1)^+
void DmdUtility::getKoopmanOperator(double input[], double output[], int rows, int cols)
{
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// [hatU1,hatSigma,hatU2]=svd(hatT(:,1:K-1),'econ')
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	int min = std::min(rows, cols - 1);

	double* leftSingVec = new double[rows * min];
	double* singVal = new double[min * min];
	double* rightSingVec = new double[min * (cols - 1)];

	double* withoutLast = new double[rows * (cols - 1)];
	DmdUtility::resizeMatrix(input, withoutLast, rows, rows, 0, cols - 2);
	// DmdUtility::printMatrix("hatT(:,2:K)", withoutLast, rows, cols - 1);

	SvdUtility::getSVD(withoutLast, rows, cols - 1, leftSingVec, singVal, rightSingVec);

	// hatT(:,2:K)
	double* withoutFirst = new double[rows * (cols - 1)];
	DmdUtility::resizeMatrix(input, withoutFirst, rows, rows, 1, cols - 1);
	// SvdUtility::printMatrix("hatT(:,2:K)", withoutFirst, rows, cols - 1);

	// hatU2	
	double* rightSingVecInv = new double[(cols - 1) * min];
	DmdUtility::transposeMatrix(rightSingVec, rightSingVecInv, min, cols - 1);
	// DmdUtility::printMatrix("hatU2", rightSingVecInv, cols - 1, min);

	// hatT(:,2:K)*hatU2
	double* xRightSingVecInv = new double[rows * min];
	DmdUtility::getMatrixMult(withoutFirst, rightSingVecInv, xRightSingVecInv, rows, cols - 1, min);
	// SvdUtility::printMatrix("xHatU2", xHatU2, n, min);

	// inv(hatSigma)
	DmdUtility::getMatrixInverse(singVal, min);
	// SvdUtility::printMatrix("hatSigma^-1", hatSigma, min, min);

	// hatT(:,2:K)*hatU2*inv(hatSigma)
	double* xSingValInv = new double[rows * min];
	DmdUtility::getMatrixMult(xRightSingVecInv, singVal, xSingValInv, rows, min, min);
	// SvdUtility::printMatrix("xHatSigmaInverse", xHatSigmaInverse, n, min);

	// hatU1'
	double* leftSingVecInv = new double[min * rows];
	DmdUtility::transposeMatrix(leftSingVec, leftSingVecInv, rows, min);
	// SvdUtility::printMatrix("hatU1^T", hatU1Transposed, min, n);

	// hatR=hatT(:,2:K)*hatU2*inv(hatSigma)*hatU1'
	DmdUtility::getMatrixMult(xSingValInv, leftSingVecInv, output, rows, min, rows);
	// DmdUtility::printMatrix("hatR", output, rows, rows);
}

// takes real matrix input (rows x cols) as double[]
// stores input^T in output (cols x rows) as double[]
void DmdUtility::transposeMatrix(double input[], double output[], int rows, int cols)
{
	for (int col = 0; col < cols; ++col)
	{
		for (int row = 0; row < rows; ++row)
		{
			output[row * cols + col] = input[col * rows + row];
		}
	}
}

// takes complex matrix input (rows x cols) as double _Complex[]
// stores input^* in output (cols x rows) as double _Complex[]
void DmdUtility::transposeMatrix(double _Complex input[], double _Complex output[], int rows, int cols)
{
	for (int col = 0; col < cols; ++col)
	{
		for (int row = 0; row < rows; ++row)
		{
			output[row * cols + col] = creal(input[col * rows + row]) - cimag(input[col * rows + row]) * csqrt(-1);
		}
	}
}

// takes square real matrix a (order x order) as double[]
// computes the inverse of a utilizing LAPACKE_dgetri
// overwrites a with a^-1
void DmdUtility::getMatrixInverse(double a[], int order)
{
	int matrix_order = LAPACK_COL_MAJOR;
	int* ipiv = new int[order];
	for (int i = 1; i <= order; ++i)
	{
		ipiv[i - 1] = i;
	}
	LAPACKE_dgetri(matrix_order, order, a, order, ipiv);
}

// takes square real matrix input (order x order) as double[]
// computes eigenvalues and eigenvectors utilizing LAPACKE_dgeev
// stores eigenvalues (order) and eigenvectors (order x order) each as double _Complex[]
void DmdUtility::getEigen(double input[], int order, double _Complex eigenvalues[], double _Complex eigenvectors[])
{
	double* eigenvaluesReal = new double[order];
	double* eigenvaluesImag = new double[order];
	double* v = new double[order * order];
	bool skip = false;

	LAPACKE_dgeev(LAPACK_COL_MAJOR, 'N', 'V', order, input, order, eigenvaluesReal, eigenvaluesImag, input, 1, v, order);

	for (int i = 0; i < order; ++i)
	{
		eigenvalues[i] = eigenvaluesReal[i] + csqrt(-1) * eigenvaluesImag[i];

		if (eigenvaluesImag[i] == 0)
		{
			for (int j = 0; j < order; ++j)
			{
				eigenvectors[j + i * order] = v[j + i * order];
			}
		}
		else if (skip)
		{
			skip = false;
		}
		else
		{
			for (int j = 0; j < order; ++j)
			{
				eigenvectors[j + i * order] = v[j + i * order] + csqrt(-1) * v[j + (i + 1) * order];
				eigenvectors[j + (i + 1) * order] = v[j + i * order] - csqrt(-1) * v[j + (i + 1) * order];
			}
			skip = true;
		}
	}
}

void DmdUtility::getDeltaOmega(double _Complex eigenvalues[], double growthRates[], double frequencies[], int size, double deltat)
{
	double _Complex* qq = new double _Complex[size];

	for (int i = 0; i < size; ++i)
	{
		qq[i] = ::clog(eigenvalues[i]);
		growthRates[i] = creal(qq[i]) / deltat;
		frequencies[i] = cimag(qq[i]) / deltat;
	}
}

void DmdUtility::getAmplitudes(double snapshots[], double _Complex eigenvalues[], double _Complex eigenvectors[], int rows, int cols, double _Complex amplitudes[])
{
	// eigenvalues=diag(MM)
	double _Complex* mm = new double _Complex[rows * rows];
	DmdUtility::setZero(mm, rows, rows);
	for (int i = 0; i < rows; ++i)
	{
		mm[i * rows + i] = eigenvalues[i];
	}

	// Mm=zeros(M*K,M)
	double _Complex* mmm = new double _Complex[rows * cols * rows];
	DmdUtility::setZero(mmm, rows * cols, rows);

	// Bb=zeros(M*K,1)
	double _Complex* bb = new double _Complex[rows * cols];
	DmdUtility::setZero(bb, rows * cols, 1);

	// MM^(k-1)
	double _Complex* mmPower = new double _Complex[rows * rows];
	DmdUtility::setZero(mmPower, rows, rows);

	// Mm(1+(k-1)*M:k*M,:)
	double _Complex* mmmMxM = new double _Complex[rows * rows];

	for (int kay = 0; kay < cols; ++kay)
	{
		DmdUtility::setZero(mmmMxM, rows, rows);
		DmdUtility::getMatrixPower(mm, mmPower, rows, kay);
		// SvdUtility::printMatrix("MM^(k-1)", mmPower, m, m);
		// Q*(MM^(k-1))
		DmdUtility::getMatrixMult(eigenvectors, mmPower, mmmMxM, rows, rows, rows);
		// SvdUtility::printMatrix("Mm(1+(k-1)*M:k*M,:)", mmmMxM, m, m);
		// Mm(1+(k-1)*M:k*M,:)=Q*(MM^(k-1))
		for (int row = 0; row < rows; ++row)
		{
			for (int col = 0; col < rows; ++col)
			{
				mmm[row + kay * rows + col * cols * rows] = mmmMxM[row + col * rows];
			}
			bb[row + kay * rows] = snapshots[row + kay * rows];
		}
	}

	double _Complex* ur = new double _Complex[rows * cols * rows];
	double* sigmar = new double[rows * rows];
	double _Complex* vrTransposed = new double _Complex[rows * rows];

	// [Ur,Sigmar,Vr]=svd(Mm,'econ')
	SvdUtility::getSVD(mmm, rows * cols, rows, ur, sigmar, vrTransposed);

	// Ur'
	double _Complex* urTransposed = new double _Complex[rows * rows * cols];
	DmdUtility::transposeMatrix(ur, urTransposed, rows * cols, rows);

	// Ur'*Bb
	double _Complex* xBb = new double _Complex[rows];
	DmdUtility::getMatrixMult(urTransposed, bb, xBb, rows, rows * cols, 1);

	// Sigmar\(Ur'*Bb)
	DmdUtility::getMatrixLeftDivision(sigmar, xBb, rows, 1);

	//Vr
	double _Complex* vr = new double _Complex[rows * rows];
	DmdUtility::transposeMatrix(vrTransposed, vr, rows, rows);

	// a=Vr*mldivide(Sigmar,Ur'*Bb)
	DmdUtility::setZero(amplitudes, rows, 1);
	DmdUtility::getMatrixMult(vr, xBb, amplitudes, rows, rows, 1);
}

// takes real matrix input (rows x cols) as double[]
// sets all entries equal zero
void DmdUtility::setZero(double input[], int rows, int cols)
{
	for (int row = 0; row < rows; ++row)
	{
		for (int col = 0; col < cols; ++col)
		{
			input[row + col * rows] = 0;
		}
	}
}

// takes complex matrix input (rows x cols) as double _Complex[]
// sets all entries equal zero
void DmdUtility::setZero(double _Complex input[], int rows, int cols)
{
	for (int row = 0; row < rows; ++row)
	{
		for (int col = 0; col < cols; ++col)
		{
			input[row + col * rows] = 0;
		}
	}
}

// takes square diagonal complex matrix input (order x order) as double _Complex[]
// stores input^exponent in output as double _Complex[]
void DmdUtility::getMatrixPower(double _Complex input[], double _Complex output[], int order, int exponent)
{
	for (int cell = 0; cell < order; ++cell)
	{
		output[cell * order + cell] = 1;
		for (int i = 0; i < exponent; ++i)
		{
			output[cell * order + cell] *= input[cell * order + cell];
			// std::cout << creal(output[cell * order + cell]) << " + " << cimag(output[cell * order + cell]) << endl;
		}
	}
}

// takes square real matrix inputA (n x n) as double[] and complex matrix input B (n x nrhs) as double _Complex[]
// computes inputA \ input B = output (n x nrhs) utilizing LAPACKE_zgesv
void DmdUtility::getMatrixLeftDivision(double inputA[], double _Complex inputB_output[], int n, int nrhs)
{
	double _Complex* inputAComplex = new double _Complex[n * n];

	for (int row = 0; row < n; ++row)
	{
		for (int col = 0; col < n; ++col)
		{
			inputAComplex[row + col * n] = inputA[row + col * n];
		}
	}

	int* ipiv = new int[n];
	for (int i = 0; i < n; ++i)
	{
		ipiv[i] = i;
	}
	LAPACKE_zgesv(LAPACK_COL_MAJOR, n, nrhs, inputAComplex, n, ipiv, inputB_output, n);
}

void DmdUtility::getDmdModes(double _Complex dmdModes[], double _Complex a[], double _Complex eigenvectors[], int rows, int cols, double leftSingVec[], double leftSingVecReduced[], double amplitudes[])
{
	// u=zeros(M:M)
	DmdUtility::setZero(dmdModes, cols, cols);

	// u(:,m)=a(m)*Q(:,m)
	for (int col = 0; col < cols; ++col)
	{
		for (int row = 0; row < cols; ++row)
		{
			dmdModes[row + col * cols] = a[col] * eigenvectors[row + col * cols];
		}
	}

	// amplitude=zeros(M,1)
	DmdUtility::setZero(amplitudes, cols, 1);

	DmdUtility::resizeMatrix(leftSingVec, leftSingVecReduced, rows, rows, 0, cols - 1);

	// aca=U*u(:,m)
	// amplitude(m)=norm(aca(:),2)/sqrt(J)
	double _Complex* aca = new double _Complex[rows];
	double _Complex* uColm = new double _Complex[cols];
	for (int col = 0; col < cols; ++col)
	{
		DmdUtility::resizeMatrix(dmdModes, uColm, cols, cols, col, col);
		DmdUtility::setZero(aca, rows, 1);
		DmdUtility::getMatrixMult(leftSingVec, uColm, aca, rows, cols, 1);
		amplitudes[col] = DmdUtility::getL2norm(aca, rows, rows) / sqrt(rows);
	}
}

int DmdUtility::getSpecComp(double _Complex dmdModes[], double growthRates[], double frequencies[], double amplitudes[], int order, double epsilon0)
{
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// UU=[u;deltas';omegas';amplitude']'
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	// transpose u
	double _Complex* dmdModesT = new double _Complex[order * order];
	DmdUtility::transposeMatrix(dmdModes, dmdModesT, order, order);
	// SvdUtility::printMatrix("uTransposed", uComplexTransposed, order, order);

	//double _Complex* deltasComplex = new double _Complex[order];
	//SvdUtility::doubleToComplex(deltas, deltasComplex, order, 1);
	// SvdUtility::printMatrix("deltas", deltasComplex, order, 1);
	//double _Complex* omegasComplex = new double _Complex[order];
	//SvdUtility::doubleToComplex(omegas, omegasComplex, order, 1);
	// SvdUtility::printMatrix("omegas", omegasComplex, order, 1);

	// UU=[u;deltas';omegas';amplitude']'
	double _Complex* uu = new double _Complex[order * (order + 3)];
	DmdUtility::concatenateVector(dmdModesT, growthRates, uu, order, order);
	DmdUtility::concatenateVector(uu, frequencies, uu, order, order + 1);
	DmdUtility::concatenateVector(uu, amplitudes, uu, order, order + 2);
	// SvdUtility::printMatrix("UU", uu, order, order + 3);

	// UU1=sortrows(UU,-(M+3))
	double _Complex* uu1 = new double _Complex[order * (order + 3)];
	DmdUtility::sortMatrix(uu, uu1, order, order + 3);
	// SvdUtility::printMatrix("UU1", uu1, order, order + 3);

	// UU=UU1'
	double _Complex* uu1Transposed = new double _Complex[(order + 3) * order];
	DmdUtility::transposeMatrix(uu1, uu1Transposed, order, order + 3);
	// SvdUtility::printMatrix("UU", uu1Transposed, order + 3, order);

	// u=UU(1:M,:)
	DmdUtility::resizeMatrix(uu1Transposed, dmdModes, order + 3, order, 0, order - 1);
	// SvdUtility::printMatrix("u", uComplex, order, order);

	// deltas=UU(M+1,:)
	// omegas=UU(M+2,:)
	// amplitude=UU(M+3,:)
	for (int col = 0; col < order; ++col)
	{
		growthRates[col] = creal(uu1Transposed[order + col * (order + 3)]);
		frequencies[col] = creal(uu1Transposed[order + 1 + col * (order + 3)]);
		amplitudes[col] = creal(uu1Transposed[order + 2 + col * (order + 3)]);
	}
	// SvdUtility::printMatrix("deltas", deltasComplex, order, 1);
	// SvdUtility::printMatrix("omegas", omegasComplex, order, 1);
	// SvdUtility::printMatrix("amplitude", amplitude, order, 1);

	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// Spectral complexity: number of DMD modes
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	int specComp = 0;
	for (int i = 0; i < order; ++i)
	{
		if (creal(amplitudes[i] / amplitudes[0]) > epsilon0)
		{
			specComp += 1;
		}
	}

	return specComp;
}

// takes complex matrix inputA (rows x cols) as double _Complex[] and real vector inputB with rows entries as double[]
// stores inputA | inputB in output (rows x cols + 1) as double _Complex[]
void DmdUtility::concatenateVector(double _Complex inputA[], double inputB[], double _Complex output[], int rows, int cols)
{
	for (int col = 0; col < cols; ++col)
	{
		for (int row = 0; row < rows; ++row)
		{
			output[row + col * rows] = inputA[row + col * rows];
		}
	}
	for (int row = 0; row < rows; ++row)
	{
		output[row + cols * rows] = inputB[row];
	}
}

// takes complex matrix input (rows x cols) as double _Complex[]
// sorts input row-wise in descending order by real part of entries of last column using selection sort
// stores sorted matrix in output as double _Complex[]
void DmdUtility::sortMatrix(double _Complex input[], double _Complex output[], int rows, int cols)
{
	double _Complex* inputCopy = new double _Complex[rows * cols];
	memcpy(inputCopy, input, rows * cols * sizeof *input);
	for (int i = 0; i < rows; ++i)
	{
		double greatestValue = creal(inputCopy[(cols - 1) * rows]);
		int greatestRow = 0;
		for (int row = 1; row < rows; ++row)
		{
			if (greatestValue < creal(inputCopy[row + (cols - 1) * rows]))
			{
				greatestValue = creal(inputCopy[row + (cols - 1) * rows]);
				greatestRow = row;
			}
		}
		for (int col = 0; col < cols; ++col)
		{
			output[i + col * rows] = inputCopy[greatestRow + col * rows];
		}
		inputCopy[greatestRow + (cols - 1) * rows] = -1;
	}
}

void DmdUtility::getDmdModesGrowthRatesFrequencies(double _Complex dmdModes[], double growthRates[], double frequencies[], double _Complex dmdModesReduced[], double growthRatesReduced[], double frequenciesReduced[], int rows, int cols, int specComp)
{
	DmdUtility::resizeMatrix(dmdModes, dmdModesReduced, rows, rows, 0, specComp - 1);

	for (int row = 0; row < specComp; ++row)
	{
		growthRatesReduced[row] = growthRates[row];
		frequenciesReduced[row] = frequencies[row];
	}
}

void DmdUtility::reconstructSnapshots(double _Complex dmdModes[], double growthRates[], double frequencies[], double leftSingVec[], double snapshotsReconst[], int rows, int cols, int specComp, int spatComp, double deltat)
{
	//hatTreconst=zeros(N,K)
	double* hatTreconst = new double[spatComp * cols];
	DmdUtility::setZero(hatTreconst, spatComp, cols);

	double* hatTreconstCol = new double[spatComp];
	for (int col = 0; col < cols; ++col)
	{
		// hatTreconst(:,k)=ContReconst_SIADS(Time(k),Time(1),u,deltas,omegas);
		DmdUtility::setZero(hatTreconstCol, spatComp, 1);
		DmdUtility::contReconst(col * deltat, 0, dmdModes, growthRates, frequencies, spatComp, specComp, hatTreconstCol);
		for (int row = 0; row < spatComp; ++row)
		{
			hatTreconst[row + col * spatComp] = hatTreconstCol[row];
		}
	}

	DmdUtility::getMatrixMult(leftSingVec, hatTreconst, snapshotsReconst, rows, spatComp, cols);
}

void DmdUtility::contReconst(double t, double t0, double _Complex dmdModes[], double growthRates[], double frequencies[], int rows, int cols, double output[])
{
	double _Complex* vv = new double _Complex[cols];
	DmdUtility::setZero(vv, cols, 1);
	for (int m = 0; m < cols; ++m)
	{
		vv[m] = cexp((growthRates[m] + csqrt(-1.0) * frequencies[m]) * (t - t0));
	}

	DmdUtility::getMatrixMult(dmdModes, vv, output, rows, cols, 1);
}

/*void DmdUtility::readCSV(string file, double output, int rows, int cols)
{

}*/

int DmdUtility::getCSVRowCount(std::string filename)
{
	int count = 0;
	std::string line;

	/* Creating input filestream */
	std::ifstream file(filename);
	while (std::getline(file, line))
		count++;

	return count;
}