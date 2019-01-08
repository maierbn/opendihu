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

  static vector<double> readCSV(string filename);

  static vector<double> readCSV(string filename, int rows);

  static void writeCSV(string filename, vector<double> values, int m, int n);

  static int getCSVRowCount(string filename);

  static int getCSVColumnCount(string filename);

};
