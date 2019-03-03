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
				cout << " ";
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
				cout << " ";
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

// takes one squre matrix input as complex double array, raises it to the exponent (>0) power and stores the resulting matrix output as complex double array
/*void SvdUtility::getMatrixPower(double _Complex input[], double _Complex output[], int order, int exponent)
{
	output = input;
	for (int i = 1; i < exponent; ++i)
	{
		cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, order, order, order, 1, output, order, input, order, 0, output, order);
	}
}*/

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

void SvdUtility::getEigen(double _Complex input[], int order, double _Complex eigenvalues[], double _Complex eigenvectors[])
{
	LAPACKE_zgeev(LAPACK_COL_MAJOR, 'N', 'V', order, input, order, eigenvalues, input, 1, eigenvectors, order);
}