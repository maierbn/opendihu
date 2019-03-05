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

// takes input data as double array and writes Sigma, U and V transposed as output parameter double arrays
// m columns, n rows
void SvdUtility::getSVD(double a[], int rows, int cols, double u[], double sigmas[], double vTransposed[], double sigma[])
{
	int min = std::min(cols, rows);
	int lda = rows, ldu = rows, ldvt = min;
	int matrix_order = LAPACK_COL_MAJOR;
	double* superb = new double[min];

	int info = LAPACKE_dgesvd(matrix_order, 's', 's', rows, cols, a, lda, sigmas, u, ldu, vTransposed, ldvt, superb);

	cout << "info: " << info << endl << endl;

	for (int row = 0; row < min; ++row)
	{
		for (int col = 0; col < min; ++col)
		{
			if (row == col)
			{
				sigma[row + col * min] = sigmas[row];
			}
			else
			{
				sigma[row + col * min] = 0;
			}

		}
	}
}

void SvdUtility::getSVD(double _Complex input[], int rows, int cols, double _Complex u[], double sigmas[], double _Complex vTransposed[], double sigma[])
{
	int min = std::min(cols, rows);
	int lda = rows, ldu = rows, ldvt = min;
	int matrix_order = LAPACK_COL_MAJOR;
	double* superb = new double[min];

	int info = LAPACKE_zgesvd(matrix_order, 's', 's', rows, cols, input, lda, sigmas, u, ldu, vTransposed, ldvt, superb);

	cout << "info: " << info << endl << endl;

	for (int row = 0; row < min; ++row)
	{
		for (int col = 0; col < min; ++col)
		{
			if (row == col)
			{
				sigma[row + col * min] = sigmas[row];
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

// returns sqrt(sum from k = size - range to size - 1 of v[k]^2)
double SvdUtility::getEuclideanNorm(double input[], int size, int range)
{
	double sum = 0;
	for (int k = size - range; k < size; ++k)
	{
		sum += input[k] * input[k];
	}
	return sqrt(sum);
}

double _Complex SvdUtility::getEuclideanNorm(double _Complex input[], int size, int range)
{
	double _Complex sum = 0;
	for (int k = size - range; k < size; ++k)
	{
		sum += creal(input[k]) * creal(input[k]) + cimag(input[k]) * cimag(input[k]);
	}
	return csqrt(sum);
}

// takes matrix a as double array, resizes it to the specified region and writes the resulting matrix in b as double array
void SvdUtility::resizeMatrix(double a[], double b[], int oldRows, int newRows, int firstCol, int lastCol)
{
	for (int col = 0; col < lastCol - firstCol + 1; ++col)
	{
		for (int row = 0; row < newRows; ++row)
		{
			b[col * newRows + row] = a[(firstCol + col) * oldRows + row];
			// cout << "column = " << col << endl << "row = " << row << endl << "a[" << (firstCol + col) * oldRows + row << "] = " << a[(firstCol + col) * oldRows + row] << endl << "b[" << col * newRows + row << "] = " << b[col * newRows + row] << endl << endl;
		}
	}
}

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

void SvdUtility::printMatrix(string name, double a[], int rows, int cols)
{
	cout << name << endl;
	for (int row = 0; row < rows; ++row)
	{
		for (int col = 0; col < cols; ++col)
		{
			if (a[col * rows + row] >= 0)
			{
				cout << " ";
			}
			cout << a[col * rows + row];
			if (col < cols - 1)
			{
				cout << "   ";
			}
		}
		cout << endl;
	}
	cout << endl;
}

void SvdUtility::printMatrix(string name, double _Complex a[], int rows, int cols)
{
	cout << name << endl;
	for (int row = 0; row < rows; ++row)
	{
		for (int col = 0; col < cols; ++col)
		{
			if (creal(a[col * rows + row]) >= 0)
			{
				cout << " ";
			}
			cout << creal(a[col * rows + row]);
			if (cimag(a[col * rows + row]) < 0)
			{
				cout << " - " << -cimag(a[col * rows + row]) << "i";
			}
			else
			{
				cout << " + " << cimag(a[col * rows + row]) << "i";
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

// takes two matrices a and b as double arrays, performs matrix multiplication and writes the resulting matrix c as double array
void SvdUtility::getMatrixProduct(double a[], double b[], double c[], int rowsA, int colsA_rowsB, int colsB, bool trans)
{
	if (trans)
	{
		cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, rowsA, colsB, colsA_rowsB, 1, a, rowsA, b, colsA_rowsB, 0, c, rowsA);
	}
	else
	{

		cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, rowsA, colsB, colsA_rowsB, 1, a, rowsA, b, colsA_rowsB, 0, c, rowsA);
	}
}

// takes two complex matrices a and b as complex double arrays, performs matrix multiplication and writes the resulting complex matrix c as complex double array
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

// takes a complex matrix input as complex double array, raises it to the eponent's power and stores the result in output
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

// takes a real matrix as double array and stores its transpose in output
void SvdUtility::transposeMatrix(double a[], double b[], int rows, int cols)
{
	for (int col = 0; col < cols; ++col)
	{
		for (int row = 0; row < rows; ++row)
		{
			b[row * cols + col] = a[col * rows + row];
			//cout << "column = " << col << endl << "row = " << row << endl << "a[" << col * rows + row << "] = " << a[col * rows + row] << endl << "b[" << row * cols + col << "] = " << b[row * cols + col] << endl << endl;
		}
	}
}

// takes a complex matrix input as complex double array and stores its transpose in output
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

// takes a real matrix a as double array, computes the inverse and stores them in a
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

// takes a complex matrix as complex double array and returns the eigenvalues and right eigenvectors as complex double arrays
void SvdUtility::getEigen(double _Complex input[], int order, double _Complex eigenvalues[], double _Complex eigenvectors[])
{
	LAPACKE_zgeev(LAPACK_COL_MAJOR, 'N', 'V', order, input, order, eigenvalues, input, 1, eigenvectors, order);
}

// takes a real matrix as double array and sets all entries equal zero
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

// takes a complex matrix as complex double array and sets all entries equal zero
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

// solves Ax = B, where x = output
void SvdUtility::getMatrixLeftDivision(double _Complex inputA[], double _Complex inputB_output[], int n, int nrhs)
{
	int* ipiv = new int[n];
	for (int i = 0; i < n; ++i)
	{
		ipiv[i] = i;
	}
	LAPACKE_zgesv(LAPACK_COL_MAJOR, n, nrhs, inputA, n, ipiv, inputB_output, n);
}

// takes a real matrix as double array and returns the same matrix as complex double array
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

// takes two matrices a (rows x cols) and b (rows x 1) and concatenates them on their row side
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

// sorts matrix rows in descending order by the positive real value of entries in last column
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
	SvdUtility::setZero(vv, 1, cols);
	for (int m = 0; m < cols; ++m)
	{
		vv[m] = cexp((deltas[m] + csqrt(-1.0) * omegas[m]) * (t - t0));
	}
	SvdUtility::getMatrixProduct(u, vv, output, rows, cols, 1);
}