#include "utility/svd_utility.h"

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
using namespace std;


void SvdUtility::getSVD(vector<double> aData)
{
  /*
   * lda = ldu = length(column)
   * ldvt= length(row)
   * s = singular values
   * superb = array to store intermediate results
   * a = matrix which we want to decomposite as 1D-array
   * 
   * */
 
   
  
  // double a[aData.size()];
  // copy(aData.begin(), aData.end(), a);
  
  int m = 6, n = 5;
  
  double a[m*n] = {
            8.79,  6.11, -9.15,  9.57, -3.49,  9.84,
            9.93,  6.91, -7.93,  1.64,  4.02,  0.15,
            9.83,  5.04,  4.86,  8.83,  9.80, -8.99,
            5.45, -0.27,  4.85,  0.74, 10.00, -6.02,
            3.16,  7.98,  3.01,  5.80,  4.27, -5.31
        };
  int lda = m, ldu = m, ldvt = n;
  int matrix_order = LAPACK_COL_MAJOR;
  //int matrix_order = LAPACK_ROW_MAJOR;
  int minmn = std::min(m,n) - 1;
  double s[n], u[ldu*m], vt[ldvt*n], superb[minmn];
	
  LAPACKE_dgesvd(matrix_order, 'a', 'a', m, n, a, lda, s, u, ldu, vt, ldvt, superb);
  
  for(int i = 0; i < ldvt*n; i++)
  {
    cout << vt[i] << endl;
  }
}

std::vector<double> SvdUtility::readCSV(string filename)
{
  ifstream data(filename);
  string line;
  vector<double> parsedCsv;
  while(getline(data, line))
  {
    stringstream lineStream(line);
    string cell;
    while(getline(lineStream,cell,','))
    {
      parsedCsv.push_back(stof(cell));
    }
  }
  return parsedCsv;
}

