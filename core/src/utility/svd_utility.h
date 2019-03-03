#pragma once

#include <Python.h>  // has to be the first included header

#include <lapacke.h>

#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <string>

using namespace std;

class SvdUtility
{

 public:
  static vector<double> getSVD(vector<double> aData, int m, int n);

  static void getSVD(double a[], int rows, int cols, double* u, double* sigmas, double* vTransposed, double* sigma);

  static vector<double> readCSV(string filename);

  static vector<double> readCSV(string filename, int rows);

  static void writeCSV(string filename, vector<double> values, int m, int n);

  static int getCSVRowCount(string filename);

  static int getCSVColumnCount(string filename);

  static double getEuclideanNorm(double input[], int size, int range);

  static void resizeMatrix(double a[], double b[], int oldRows, int newRows, int firstCol, int lastCol);

  static void printMatrix(string name, double a[], int rows, int columns);

  static void printMatrix(string name, double _Complex a[], int rows, int columns);

  static void getMatrixProduct(double a[], double b[], double c[], int rowsA, int colsA_rowsB, int colsB, bool trans);

  //static void getMatrixPower(double _Complex input[], double _Complex output[], int order, int exponent);

  static void transposeMatrix(double a[], double b[], int rows, int cols);

  static void getMatrixInverse(double a[], int order);

  static void getEigen(double _Complex input[], int order, double _Complex eigenvalues[], double _Complex eigenvectors[]);
};
