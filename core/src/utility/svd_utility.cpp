#include "utility/svd_utility.h"

#include <cblas.h>

#include <iostream>
#include <vector>
#include <string>
#include <fstream>

using namespace std;

// applies SVD and returns V transposed
// m columns and n rows
std::vector<double> SvdUtility::getSVD(vector<double> aData, int m, int n)
{
	/*
	 * lda = ldu = length(column)
	 * ldvt= length(row)
	 * s = singular values
	 * superb = array to store intermediate results
	 * a = matrix which we want to decomposite as 1D-array
	 *
	 * */



	double* a = new double[aData.size()];
	copy(aData.begin(), aData.end(), a);
	// Spalten    Zeilen
	// int m = 6, n = 5;

	/* double a[m*n] = {
			  8.79,  6.11, -9.15,  9.57, -3.49,  9.84,
			  9.93,  6.91, -7.93,  1.64,  4.02,  0.15,
			  9.83,  5.04,  4.86,  8.83,  9.80, -8.99,
			  5.45, -0.27,  4.85,  0.74, 10.00, -6.02,
			  3.16,  7.98,  3.01,  5.80,  4.27, -5.31
		  };
	**/
	int lda = m, ldu = m, ldvt = n;
	int matrix_order = LAPACK_COL_MAJOR;
	//int matrix_order = LAPACK_ROW_MAJOR;
	int minmn = std::min(m, n) - 1;
	double* s = new double[n];
	double* u = new double[ldu*m];
	double* vt = new double[ldvt*n];
	double* superb = new double[minmn];

	int info = LAPACKE_dgesvd(matrix_order, 'a', 'a', m, n, a, lda, s, u, ldu, vt, ldvt, superb);

	cout << "info: " << info << endl;

	for (int i = 0; i < ldvt*n; i++)
	{
		cout << vt[i] << endl;
	}

	return std::vector<double>(vt, vt + sizeof vt / sizeof vt[0]);;
}

// takes real matrix input (rows x cols) as double[] in column major order
// performs singular-value decomposition on input utilizing LAPACKE_dgesvd
// stores the left-singular vectors (column-wise leftSingularVectors), singular values (singularValues as vector, diagonal entries of sigma as matrix) and right-singular vectors (row-wise rightSingularVectorsTransposed) as double arrays
void SvdUtility::getSVD(double input[], int rows, int cols, double leftSingularVectors[], double singularValues[], double rightSingularVectorsTransposed[], double sigma[])
{
	int min = std::min(cols, rows);
	double* superb = new double[min];

	int info = LAPACKE_dgesvd(LAPACK_COL_MAJOR, 's', 's', rows, cols, input, rows, singularValues, leftSingularVectors, rows, rightSingularVectorsTransposed, min, superb);

	cout << "info: " << info << endl << endl;

	// build diagonal matrix sigma from vector singularValues
	for (int row = 0; row < min; ++row)
	{
		for (int col = 0; col < min; ++col)
		{
			if (row == col)
			{
				sigma[row + col * min] = singularValues[row];
			}
			else
			{
				sigma[row + col * min] = 0;
			}

		}
	}
}

// takes complex matrix input (rows x cols) as double _Complex[] in column major order
// performs singular-value decomposition on input utilizing LAPACKE_zgesvd
// stores the left-singular vectors (column-wise leftSingularVectors), singular values (singularValues as vector, diagonal entries of sigma as matrix) and right-singular vectors (row-wise rightSingularVectorsTransposed) as double _Complex arrays
void SvdUtility::getSVD(double _Complex input[], int rows, int cols, double _Complex leftSingularVectors[], double singularValues[], double _Complex rightSingularVectorsTransposed[], double sigma[])
{
	int min = std::min(cols, rows);
	double* superb = new double[min];

	int info = LAPACKE_zgesvd(LAPACK_COL_MAJOR, 's', 's', rows, cols, input, rows, singularValues, leftSingularVectors, rows, rightSingularVectorsTransposed, min, superb);

	cout << "info: " << info << endl << endl;

	// build diagonal matrix sigma from vector singularValues
	for (int row = 0; row < min; ++row)
	{
		for (int col = 0; col < min; ++col)
		{
			if (row == col)
			{
				sigma[row + col * min] = singularValues[row];
			}
			else
			{
				sigma[row + col * min] = 0;
			}

		}
	}
}

// reads CSV cell by cell as vector
std::vector<double> SvdUtility::readCSV(string filename)
{
	ifstream data(filename);
	string line;
	vector<double> parsedCsv;
	int i, j = 0;
	while (getline(data, line))
	{
		stringstream lineStream(line);
		i++;
		string cell;
		j = 0;
		while (getline(lineStream, cell, ','))
		{
			parsedCsv.push_back(stof(cell));
			j++;
		}
		// cout << i << endl;
		// cout << j << endl;
	}
	return parsedCsv;
}

// reads specified rows of CSV cell by cell as vector
std::vector<double> SvdUtility::readCSV(string filename, int rows)
{
	ifstream data(filename);
	string line;
	vector<double> parsedCsv;
	for (int i = 0; i < rows; i++)
	{
		getline(data, line);
		stringstream lineStream(line);
		string cell;
		while (getline(lineStream, cell, ','))
		{
			parsedCsv.push_back(stof(cell));
		}
	}
	return parsedCsv;
}

// writes vector cell by cell as CSV
void SvdUtility::writeCSV(string filename, std::vector<double> values, int m, int n)
{
	ofstream data;
	data.open(filename);
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			data << std::to_string(values[i*m + j]) << ",";
		}
		data << "\n";
	}
	data.close();
}

// returns number of rows of CSV
int SvdUtility::getCSVRowCount(string filename)
{
	ifstream data(filename);
	string line;
	int i = 0;
	while (getline(data, line))
	{
		stringstream lineStream(line);
		i++;
	}
	return i;
}

// returns number of columns of CSV
int SvdUtility::getCSVColumnCount(string filename)
{
	ifstream data(filename);
	string line;
	vector<double> parsedCsv;
	int j = 0;
	if (getline(data, line))
	{
		stringstream lineStream(line);
		string cell;
		while (getline(lineStream, cell, ','))
		{
			parsedCsv.push_back(stof(cell));
			j++;
		}
	}
	return j;
}

// takes real vector input with size entries as double[]
// returns the Euclidean norm of vector input over the specified range of entries as double
// sqrt(sum from k = size - range to size - 1 of v[k]^2)
double SvdUtility::getEuclideanNorm(double input[], int size, int range)
{
	double sum = 0;
	for (int k = size - range; k < size; ++k)
	{
		sum += input[k] * input[k];
	}
	return sqrt(sum);
}

// takes complex vector input with size entries as double _Complex[]
// returns the Euclidean norm of vector input over the specified range of entries as double _Complex
// sqrt(sum from k = size - range to size - 1 of Re(v[k])^2 + Im(v[k])^2)
double _Complex SvdUtility::getEuclideanNorm(double _Complex input[], int size, int range)
{
	double _Complex sum = 0;
	for (int k = size - range; k < size; ++k)
	{
		sum += creal(input[k]) * creal(input[k]) + cimag(input[k]) * cimag(input[k]);
	}
	return csqrt(sum);
}

// takes real matrix input with oldRows rows as double[]
// stores the entries of input in real matrix output with newRows rows and lastCol - firstCol + 1 columns as double[] starting with input[0, firstCol] and ending with input[newRows - 1, lastCol]
void SvdUtility::resizeMatrix(double input[], double output[], int oldRows, int newRows, int firstCol, int lastCol)
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
void SvdUtility::resizeMatrix(double _Complex input[], double _Complex output[], int oldRows, int newRows, int firstCol, int lastCol)
{
	for (int col = 0; col < lastCol - firstCol + 1; ++col)
	{
		for (int row = 0; row < newRows; ++row)
		{
			output[col * newRows + row] = input[(firstCol + col) * oldRows + row];
		}
	}
}

// takes real matrix input (rows x cols) as double[]
// prints name and input entry-wise
void SvdUtility::printMatrix(string name, double input[], int rows, int cols)
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
void SvdUtility::printMatrix(string name, double _Complex input[], int rows, int cols)
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

// takes two real matrices inputA (rowsA x colsA_rowsB) and inputB (colsA_rowsB x colsB) each as double[]
// performs matrix multiplication utilizing cblas_dgemm
// stores if (trans) {inputA * inputB^T} else {inputA * inputB} in output (rowsA x colsB) as double[]
void SvdUtility::getMatrixProduct(double inputA[], double inputB[], double output[], int rowsA, int colsA_rowsB, int colsB, bool trans)
{
	if (trans)
	{
		cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, rowsA, colsB, colsA_rowsB, 1, inputA, rowsA, inputB, colsA_rowsB, 0, output, rowsA);
	}
	else
	{

		cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, rowsA, colsB, colsA_rowsB, 1, inputA, rowsA, inputB, colsA_rowsB, 0, output, rowsA);
	}
}

// takes two complex matrices inputA (rowsA x colsA_rowsB) and inputB (colsA_rowsB x colsB) each as double _Complex[]
// stores inputA * inputB in output (rowsA x colsB) as double _Complex[]
void SvdUtility::getMatrixProduct(double _Complex inputA[], double _Complex inputB[], double _Complex output[], int rowsA, int colsA_rowsB, int colsB)
{
	// cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, rowsA, colsB, colsA_rowsB, 1, inputA, rowsA, inputB, colsA_rowsB, 0, output, rowsA);

	for (int col = 0; col < colsB; ++col)
	{
		for (int row = 0; row < rowsA; ++row)
		{
			for (int cell = 0; cell < colsA_rowsB; ++cell)
			{
				output[row + col * rowsA] += inputA[row + cell * rowsA] * inputB[col * colsA_rowsB + cell];
			}
		}
	}
}

// takes square diagonal complex matrix input (order x order) as double _Complex[]
// stores input^exponent in output as double _Complex[]
void SvdUtility::getMatrixPower(double _Complex input[], double _Complex output[], int order, int exponent)
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

// takes real matrix input (rows x cols) as double[]
// stores input^T in output (cols x rows) as double[]
void SvdUtility::transposeMatrix(double input[], double output[], int rows, int cols)
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
// stores input^T in output (cols x rows) as double _Complex[]
void SvdUtility::transposeMatrix(double _Complex input[], double _Complex output[], int rows, int cols)
{
	for (int col = 0; col < cols; ++col)
	{
		for (int row = 0; row < rows; ++row)
		{
			output[row * cols + col] = input[col * rows + row];
		}
	}
}

// takes square real matrix a (order x order) as double[]
// computes the inverse of a utilizing LAPACKE_dgetri
// overwrites a with a^-1
void SvdUtility::getMatrixInverse(double a[], int order)
{
	int matrix_order = LAPACK_COL_MAJOR;
	int* ipiv = new int[order];
	for (int i = 1; i <= order; ++i)
	{
		ipiv[i - 1] = i;
	}
	LAPACKE_dgetri(matrix_order, order, a, order, ipiv);
}

// takes square complex matrix input (order x order) as double _Complex[]
// computes eigenvalues and eigenvectors utilizing LAPACKE_zgeev
// stores eigenvalues (order) and eigenvectors (order x order) each as double _Complex[] array
void SvdUtility::getEigen(double _Complex input[], int order, double _Complex eigenvalues[], double _Complex eigenvectors[])
{
	LAPACKE_zgeev(LAPACK_COL_MAJOR, 'N', 'V', order, input, order, eigenvalues, input, 1, eigenvectors, order);
}

// takes real matrix input (rows x cols) as double[]
// sets all entries equal zero
void SvdUtility::setZero(double input[], int rows, int cols)
{
	for (int row = 0; row < rows; ++row)
	{
		for (int col = 0; col < cols; ++col)
		{
			input[row + col * rows] = 0.0;
		}
	}
}

// takes complex matrix input (rows x cols) as double _Complex[]
// sets all entries equal zero
void SvdUtility::setZero(double _Complex input[], int rows, int cols)
{
	for (int row = 0; row < rows; ++row)
	{
		for (int col = 0; col < cols; ++col)
		{
			input[row + col * rows] = 0.0;
		}
	}
}

// takes square complex matrix inputA (n x n) and complex matrix input B (n x nrhs) each as double _Complex[]
// computes inputA \ input B = output (n x nrhs) utilizing LAPACKE_zgesv
void SvdUtility::getMatrixLeftDivision(double _Complex inputA[], double _Complex inputB_output[], int n, int nrhs)
{
	int* ipiv = new int[n];
	for (int i = 0; i < n; ++i)
	{
		ipiv[i] = i;
	}
	LAPACKE_zgesv(LAPACK_COL_MAJOR, n, nrhs, inputA, n, ipiv, inputB_output, n);
}

// takes real matrix input (rows x cols) as double[]
// stores each entry of input as real value in complex matrix output (rows x cols) as double _Complex[]
void SvdUtility::doubleToComplex(double input[], double _Complex output[], int rows, int cols)
{
	for (int row = 0; row < rows; ++row)
	{
		for (int col = 0; col < cols; ++col)
		{
			output[row + col * rows] = input[row + col * rows];
		}
	}
}

// takes complex matrix inputA (rows x cols) and complex vector inputB with rows entries each as double _Complex[]
// stores inputA | inputB in output (rows x cols + 1) as double _Complex[]
void SvdUtility::concatenateMatrices(double _Complex inputA[], double _Complex inputB[], double _Complex output[], int rows, int cols)
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
void SvdUtility::sortMatrix(double _Complex input[], double _Complex output[], int rows, int cols)
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

void SvdUtility::contReconst(double t, double t0, double _Complex u[], double _Complex deltas[], double _Complex omegas[], int rows, int cols, double _Complex output[])
{
	double _Complex* vv = new double _Complex[cols];
	SvdUtility::setZero(vv, cols, 1);
	for (int m = 0; m < cols; ++m)
	{
		vv[m] = cexp((deltas[m] + csqrt(-1.0) * omegas[m]) * (t - t0));
	}
	SvdUtility::getMatrixProduct(u, vv, output, rows, cols, 1);
}