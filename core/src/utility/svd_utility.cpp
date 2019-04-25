#include "utility/svd_utility.h"

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
using namespace std;


std::vector<double> SvdUtility::getSVD(vector<double> &aData, int m, int n)
{
  /*
   * lda = ldu = length(column)
   * ldvt= length(row)
   * s = singular values
   * superb = array to store intermediate results
   * a = matrix which we want to decomposite as 1D-array
   * 
   * */
 
   
  
  //double a[aData.size()];
  //copy(aData.begin(), aData.end(), a);
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
  int minmn = std::min(m,n) - 1;
  std::vector<double> s(n), u(ldu*m), vt(ldvt*n), superb(minmn);
	
  LAPACKE_dgesvd(matrix_order, 'a', 'a', m, n, aData.data(), lda, s.data(), u.data(), ldu, vt.data(), ldvt, superb.data());
  
  for(int i = 0; i < ldvt*n; i++)
  {
    cout << vt[i] << endl;
  }
  
  //return std::vector<double>(vt, vt + sizeof vt / sizeof vt[0]);  // <-- ??
  return vt;
}

std::vector<double> SvdUtility::readCSV(string filename)
{
  ifstream data(filename);
  string line;
  vector<double> parsedCsv;
  int i, j = 0;
  while(getline(data, line))
  {
    stringstream lineStream(line);
    i++;
    string cell;
    j = 0;
    while(getline(lineStream,cell,','))
    {
      parsedCsv.push_back(stof(cell));
      j++;
    }
    // cout << i << endl;
    // cout << j << endl;
  }
  return parsedCsv;
}

std::vector<double> SvdUtility::readCSV(string filename, int rows)
{
  ifstream data(filename);
  string line;
  vector<double> parsedCsv;
  for(int i = 0; i < rows; i++)
  {
    getline(data, line);
    stringstream lineStream(line);
    string cell;
    while(getline(lineStream,cell,','))
    {
      parsedCsv.push_back(stof(cell));
    }
  }
  return parsedCsv;
}

void SvdUtility::writeCSV(string filename, std::vector<double> &values, int m, int n)
{
   ofstream data;
   data.open(filename);
   for(int i = 0; i < m; i++)
   {
     for(int j = 0; j < n; j++)
     {
       data << std::to_string(values[i*m + j]) << ",";
     }
     data << "\n";
   }
   data.close();
}

int SvdUtility::getCSVRowCount(string filename)
{
  ifstream data(filename);
  string line;
  int i = 0;
  while(getline(data, line))
  {
    stringstream lineStream(line);
    i++;
  }
  return i;
}

int SvdUtility::getCSVColumnCount(string filename)
{
  ifstream data(filename);
  string line;
  vector<double> parsedCsv;
  int j = 0;
  if(getline(data, line))
  {
    stringstream lineStream(line);
    string cell;
    while(getline(lineStream,cell,','))
    {
      parsedCsv.push_back(stof(cell));
      j++;
    }
  }
  return j;
}
