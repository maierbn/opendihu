#include <iostream>
#include <cstdlib>
#include <vector>
#include "opendihu.h"
#include "utility/svd_utility.h"

int main(int argc, char *argv[])
{
  std::vector<double> parsedCSV = SvdUtility::readCSV("./out/data.csv");
  for(std::vector<double>::iterator it = parsedCSV.begin(); it!=parsedCSV.end(); ++it)
  {    
    std::cout << ' ' << *it;
  }
  
  std::cout << std::endl;

  int rowCount = SvdUtility::getCSVRowCount("./out/data.csv");
  int columnCount = SvdUtility::getCSVColumnCount("./out/data.csv");
  std::cout << "rowCount: " << rowCount << " columnCount: " << columnCount << std::endl;
  std::cout << "getSVD: " << std::endl;
  std::vector<double> result = SvdUtility::getSVD(parsedCSV, columnCount, rowCount);
  
  for(std::vector<double>::iterator it = result.begin(); it!=result.end(); ++it)
  {    
    std::cout << ' ' << *it;
  }
  
  std::cout << "/n" << std::endl;
  
  std::cout << "size Vt: " << result.size() << std::endl;
  SvdUtility::writeCSV("./out/SVDresult.csv", result, 51, 51);
  return EXIT_SUCCESS;
}
