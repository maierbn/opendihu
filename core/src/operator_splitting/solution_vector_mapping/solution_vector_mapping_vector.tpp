#include "operator_splitting/solution_vector_mapping/solution_vector_mapping_vector.h"

#include <vector>
#include <tuple>
#include "easylogging++.h"

/** Transfer between two vectors of any type
 */
template<typename TransferableSolutionDataType1, typename TransferableSolutionDataType2>
void SolutionVectorMapping<
  std::vector<TransferableSolutionDataType1>,
  std::vector<TransferableSolutionDataType2>
>::transfer(const std::vector<TransferableSolutionDataType1> &transferableSolutionData1,
            const std::vector<TransferableSolutionDataType2> &transferableSolutionData2)
{
  int nTransferableVariables = std::min(transferableSolutionData1.size(), transferableSolutionData2.size());
  if (transferableSolutionData1.size() != transferableSolutionData2.size())
  {
    LOG(ERROR) << "Trying to transfer data from " << transferableSolutionData1.size() << " variables to " << transferableSolutionData2.size() << ", number has to be equal. "
      << "Now only using the first " << nTransferableVariables << " variable" << (nTransferableVariables != 1? "s" : "") << ". Types: " << std::endl
      << typeid(TransferableSolutionDataType1).name() << std::endl << typeid(TransferableSolutionDataType2).name();
  }

  for (int i = 0; i < nTransferableVariables; i++)
  {
    SolutionVectorMapping<TransferableSolutionDataType1,TransferableSolutionDataType2>::transfer(
      transferableSolutionData1[i],
      transferableSolutionData2[i]
    );
  }
}
