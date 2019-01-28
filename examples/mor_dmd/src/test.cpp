#include <iostream>
#include <cstdlib>
#include <vector>
#include "opendihu.h"
#include "utility/svd_utility.h"

int main(int argc, char *argv[])
{
  // input data
  //std::vector<double> parsedCsv = SvdUtility::readCSV("./out/data.csv");
  //double v[parsedCsv.size()];
  //copy(parsedCsv.begin(), parsedCsv.end(), v);

  //int rowCount = SvdUtility::getCSVRowCount("./out/data.csv");
  //int columnCount = SvdUtility::getCSVColumnCount("./out/data.csv");
  
  int rowCount = 5;
  int columnCount = 6;

  std::cout << "rowCount: " << rowCount << ", columnCount: " << columnCount << endl << endl;

  double v[rowCount*columnCount] = {
            8.79,  6.11, -9.15,  9.57, -3.49,  9.84,
            9.93,  6.91, -7.93,  1.64,  4.02,  0.15,
            9.83,  5.04,  4.86,  8.83,  9.80, -8.99,
            5.45, -0.27,  4.85,  0.74, 10.00, -6.02,
            3.16,  7.98,  3.01,  5.80,  4.27, -5.31
        };

  // n=length(sigmas);
  int n = std::min(rowCount, columnCount);


  // [U,Sigma,T]=svd(V,'econ');
  double u[columnCount*columnCount];
  double sigma[n];
  double t[rowCount*rowCount];
  SvdUtility::getSVD(v, columnCount, rowCount, u, sigma, t);

  std::cout << endl;
  std::cout << "Singular values" << endl;
  for(int i=0;i<n;++i)
  {
    std::cout << sigma[i] << " ";
  }
  std::cout << endl << endl;

  std::cout << "Left singular values (stored columnwise)" << endl;
  for(int i=0;i<columnCount;++i)
  {
    for(int j=0;j<columnCount;++j)
    {
      std::cout << u[i*columnCount+j] << " ";
    }
    std::cout << endl;
  }
  std::cout << endl;

  std::cout << "Right singular values (stored rowwise)" << endl;
  for(int i=0;i<rowCount;++i)
  {
    for(int j=0;j<rowCount;++j)
    {
      std::cout << t[i*rowCount+j] << " ";
    }
    std::cout << endl;
  }
  std::cout << endl;

  std::cout << "SVD done" << endl;

  // sigmas=diag(Sigma);
  // NormS=norm(sigmas,2);
  //double sigmas [n];
  double normS = 0;
  /*for(int i = 0; i < n; ++i)
  {
    sigmas[i] = sigma[i + columnCount*i];
    std::cout << sigmas[i];
    normS += pow(sigmas[i], 2.0);
  }*/
  normS = pow(normS, 0.5);

  std::cout << normS << endl;

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
