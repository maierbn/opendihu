#include <iostream>
#include <cstdlib>
#include <vector>
#include "opendihu.h"
#include "utility/svd_utility.h"

int main(int argc, char *argv[])
{
  // input data
  std::vector<double> parsedCsv = SvdUtility::readCSV("./out/data.csv");
  double v[parsedCsv.size()];
  copy(parsedCsv.begin(), parsedCsv.end(), v);

  int rowCount = SvdUtility::getCSVRowCount("./out/data.csv");
  int columnCount = SvdUtility::getCSVColumnCount("./out/data.csv");
  std::cout << "rowCount: " << rowCount << " columnCount: " << columnCount << std::endl;

  // [U,Sigma,T]=svd(V,'econ');
  double u [rowCount*rowCount];
  double sigma [rowCount*columnCount];
  double t [columnCount*columnCount];
  SvdUtility::getSVD(v, columnCount, rowCount, u, sigma, t);

  // n=length(sigmas);
  int n = std::min(rowCount, columnCount);

  // sigmas=diag(Sigma);
  // NormS=norm(sigmas,2);
  double sigmas [n];
  double normS = 0;
  for(int i = 0; i < n; ++i)
  {
    sigmas[i] = sigma[i + columnCount*i];
    std::cout << sigmas[i];
    normS += pow(sigmas[i], 2.0);
  }
  normS = pow(normS, 0.5);

  std::cout << normS;

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
