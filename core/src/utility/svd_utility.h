#pragma once

#include <Python.h>  // has to be the first included header

#ifdef HAVE_LAPACK
#include <lapacke.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <string>

using namespace std;

class SvdUtility
{

 public:
  
  static void getSVD(double input[], int rows, int cols, double leftSingVec[]);

  static void getSVD(double input[], int rows, int cols, double leftSingVec[], double sigma[], double rightSingVecT[]);

  static void getSVD(double _Complex input[], int rows, int cols, double _Complex leftSingVec[], double sigma[], double _Complex rightSingVecT[]);

  static void reconstructSnapshots(int rows, int cols, double leftSingVec[], double sigma[], double rightSingVecT[], double output[]);
  
  static void printMatrix(std::string name, double input[], int rows, int cols);
  
  static vector<double> readCSV(string filename);

  static vector<double> readCSV(string filename, int rows);

  static void writeCSV(string filename, vector<double> values, int m, int n);
  
  static void writeCSV(string filename, double values[], int m, int n, bool columnWise = false);
  
  static int getCSVRowCount(string filename);

  static int getCSVColumnCount(string filename);
};


