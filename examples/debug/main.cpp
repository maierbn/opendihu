#include <iostream>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <array>

#include "easylogging++.h"
#include "semt/Semt.h"

#include "opendihu.h"

int main(int argc, char *argv[])
{
  std::string pythonConfig = R"(
config = {
  "nElements": [2,2],
}
)";

  DihuContext settings(argc, argv, pythonConfig);
  
  typedef FunctionSpace::FunctionSpace<
    Mesh::StructuredRegularFixedOfDimension<1>,
    BasisFunction::LagrangeOfOrder<1>
  > FunctionSpaceType;
  
  
  std::shared_ptr<FunctionSpaceType> functionSpace = settings.meshManager()->createFunctionSpace<FunctionSpaceType>("test", settings.getPythonConfig());
  
  const int nComponentsTest = 1;

  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nComponentsTest>> test = functionSpace->template createFieldVariable<nComponentsTest>("test");
  
  //Vec globalVector = test->valuesGlobal();
  //Vec localVector = test->valuesLocal();
  
  //test->zeroEntries();
  //test->zeroGhostBuffer();
  test->setRepresentationGlobal();
  test->startGhostManipulation();
  LOG(DEBUG) << "inititalized, test: " << *test;


  int ownRankNo = functionSpace->meshPartition()->rankSubset()->ownRankNo();
  if (ownRankNo == 0)
  {
    std::array<int,2> indices({0,1});
    std::array<double,2> values({0.1,1.1});
      
    // here VecSetValues is used, not VecSetValuesLocal, because it is a plain vector
    test->setValues<2>(indices, values, INSERT_VALUES);
  }
  else
  {
    
    std::array<int,2> indices({0,1});
    std::array<double,2> values({0.2,1.2});
          
    test->setValues<2>(indices, values, INSERT_VALUES);
  }

  LOG(DEBUG) << "before finishGhostManipulation, test: " << *test;
  test->finishGhostManipulation();


  LOG(DEBUG) << "after finishGhostManipulation, test: " << *test;


  test->debug();
  
  return EXIT_SUCCESS;
}
