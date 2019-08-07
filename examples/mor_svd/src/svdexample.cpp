#include <iostream>
#include <cstdlib>
#include <vector>
#include "opendihu.h"
#include "utility/svd_utility.h"
#include "utility/dmd_utility.h"

int main(int argc, char *argv[])
{
  std::vector<double> parsedCSV = SvdUtility::readCSV("./out/snapshots.csv");
  for(std::vector<double>::iterator it = parsedCSV.begin(); it!=parsedCSV.end(); ++it)
  {    
    std::cout << ' ' << *it;
  } 
  std::cout << std::endl;
  
  // parsedCSV is the transpose of the snapshots matrix
  int rowsSnapshots = SvdUtility::getCSVColumnCount("./out/snapshots.csv");
  int columnsSnapshots = SvdUtility::getCSVRowCount("./out/snapshots.csv");
  std::cout << "rowsSnapshots: " << rowsSnapshots << " columnsSnapshots: " << columnsSnapshots << std::endl;
  
  double* snapshots = new double[rowsSnapshots * columnsSnapshots];
  std::copy(parsedCSV.begin(), parsedCSV.end(), snapshots);  
  int rankSnapshots = std::min(rowsSnapshots , columnsSnapshots);
  
  // make a copy before the values of the input are changed by the lapack SVD
  double* input=new double[rowsSnapshots * columnsSnapshots]; 
  std::copy(snapshots,snapshots+rowsSnapshots * columnsSnapshots,input); 
  
  double* leftSingVec = new double[rowsSnapshots * rankSnapshots];
  double* sigma = new double[rankSnapshots * rankSnapshots];
  double* rightSingVecT = new double[rankSnapshots*columnsSnapshots];
  
  std::cout << "getSVD: " << std::endl;
  //SvdUtility::getSVD(snapshots, rowsSnapshots, columnsSnapshots, leftSingVec);
  SvdUtility::getSVD(input, rowsSnapshots, columnsSnapshots, leftSingVec, sigma, rightSingVecT);
  
  double* snapshots_reconst = new double[rowsSnapshots * columnsSnapshots];
  SvdUtility::reconstructSnapshots(rowsSnapshots, columnsSnapshots, leftSingVec, sigma, rightSingVecT,snapshots_reconst);
  SvdUtility::writeCSV("./out/snapshots_svdReconst.csv", snapshots_reconst, columnsSnapshots, rowsSnapshots); // transpose of the snapshots are reconstructed
  SvdUtility::writeCSV("./out/SVDresult.csv", leftSingVec, rowsSnapshots, rankSnapshots );
  
  double* diff = new double[rowsSnapshots * columnsSnapshots];
  for (int row = 0; row < rowsSnapshots; ++row)
  {
    for (int col = 0; col < columnsSnapshots; ++col)
    {
      diff[row * columnsSnapshots + col] = snapshots[row * columnsSnapshots + col] - snapshots_reconst[row * columnsSnapshots + col];
      //std::cout << "snapshot: " << snapshots[row * columnsSnapshots + col] << ", recons: " << snapshots_reconst[row * columnsSnapshots + col] << ", diff: " << diff[row * columnsSnapshots + col] << std::endl;      
    }
  }
  double normSnapshots = DmdUtility::getNorm(snapshots, rowsSnapshots, columnsSnapshots);
  double rms = DmdUtility::getNorm(diff, rowsSnapshots, columnsSnapshots) / normSnapshots;
  std::cout << "RelativeerrorRMS = " << rms << endl;
  
  
  return EXIT_SUCCESS;
}
