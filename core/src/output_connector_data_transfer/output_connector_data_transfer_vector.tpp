#include "output_connector_data_transfer/output_connector_data_transfer_vector.h"

#include <vector>
#include <tuple>
#include "easylogging++.h"

/** Transfer between two vectors of any type
 */
template<typename OutputConnectorDataType1, typename OutputConnectorDataType2>
void SolutionVectorMapping<
  std::vector<std::shared_ptr<OutputConnectorDataType1>>,
  std::vector<std::shared_ptr<OutputConnectorDataType2>>
>::transfer(const std::shared_ptr<std::vector<std::shared_ptr<OutputConnectorDataType1>>> transferableSolutionData1,
            std::shared_ptr<std::vector<std::shared_ptr<OutputConnectorDataType2>>> transferableSolutionData2,
            OutputConnection &outputConnection)
{
  LOG(DEBUG) << "transfer vector (1)";

  VLOG(1) << "Solution vector mapping (output_connector_data_transfer_vector.tpp)";
  int nTransferableVariables = std::min(transferableSolutionData1->size(), transferableSolutionData2->size());
  if (transferableSolutionData1->size() != transferableSolutionData2->size())
  {
    LOG(ERROR) << "Trying to transfer data from " << transferableSolutionData1->size() << " variables to " << transferableSolutionData2->size() << ", number has to be equal. "
      << "Now only using the first " << nTransferableVariables << " variable" << (nTransferableVariables != 1? "s" : "") << ". Types: " << std::endl
      << StringUtility::demangle(typeid(OutputConnectorDataType1).name()) << std::endl
      << StringUtility::demangle(typeid(OutputConnectorDataType2).name());
  }

  // loop over the two vectors A and B and transfer A[i] -> B[i] for every component
  for (int i = 0; i < nTransferableVariables; i++)
  {
    SolutionVectorMapping<OutputConnectorDataType1,OutputConnectorDataType2>::transfer(
      (*transferableSolutionData1)[i],
      (*transferableSolutionData2)[i],
      outputConnection
    );
  }

  VLOG(1) << "at the end of std::vector transfer: transferableSolutionData2: " << (*transferableSolutionData2);
}

/** Transfer between a normal entry and a vector, the first entry of the vector is used
 */
template<typename OutputConnectorDataType1, typename OutputConnectorDataType2>
void SolutionVectorMapping<
  OutputConnectorDataType1,
  std::vector<std::shared_ptr<OutputConnectorDataType2>>
>::transfer(const std::shared_ptr<OutputConnectorDataType1> transferableSolutionData1,
            std::shared_ptr<std::vector<std::shared_ptr<OutputConnectorDataType2>>> transferableSolutionData2,
            OutputConnection &outputConnection)
{
  LOG(DEBUG) << "transfer vector (2)";
  VLOG(1) << "Solution vector mapping (output_connector_data_transfer_vector.tpp)";

  if (transferableSolutionData2->empty())
  {
    LOG(ERROR) << "Trying to transfer data from a single entry to an empty vector. Types: " << std::endl
      << StringUtility::demangle(typeid(OutputConnectorDataType1).name()) << std::endl
      << StringUtility::demangle(typeid(OutputConnectorDataType2).name());
    return;
  }

  SolutionVectorMapping<OutputConnectorDataType1,OutputConnectorDataType2>::transfer(
    transferableSolutionData1,
    (*transferableSolutionData2)[0],
    outputConnection
  );
}

/** Transfer between a normal entry and a vector, the first entry of the vector is used
 */
template<typename OutputConnectorDataType1, typename OutputConnectorDataType2>
void SolutionVectorMapping<
  std::vector<std::shared_ptr<OutputConnectorDataType1>>,
  OutputConnectorDataType2
>::transfer(const std::shared_ptr<std::vector<std::shared_ptr<OutputConnectorDataType1>>> transferableSolutionData1,
            std::shared_ptr<OutputConnectorDataType2> transferableSolutionData2,
            OutputConnection &outputConnection)
{
  LOG(DEBUG) << "transfer vector (3)";
  VLOG(1) << "Solution vector mapping (output_connector_data_transfer_vector.tpp)";

  if (transferableSolutionData1->empty())
  {
    LOG(ERROR) << "Trying to transfer data from a single entry to an empty vector. Types: " << std::endl
      << StringUtility::demangle(typeid(OutputConnectorDataType1).name()) << std::endl
      << StringUtility::demangle(typeid(OutputConnectorDataType2).name());
    return;
  }

  SolutionVectorMapping<OutputConnectorDataType1,OutputConnectorDataType2>::transfer(
    (*transferableSolutionData1)[0],
    transferableSolutionData2,
    outputConnection
  );
}
