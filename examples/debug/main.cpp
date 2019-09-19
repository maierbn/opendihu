
#include <Vc/Vc>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <array>
#include <cmath>

#include "easylogging++.h"
#include "semt/Semt.h"

#include "opendihu.h"

using namespace std;

struct FiberPointData
{
  Vc::double_v values[8];
};
Vc_DECLARE_ALLOCATOR(FiberPointData)

int f(int a)
{
  double v = 0;
  for(long long int j = 0; j < 1e8; j++)
  {
    if (j % 17 == 5 || j % 13 == 2 || j%7 == 6)
      v += 2*a;
    if (j%2 == 0)
      v *= 2.;
    if (j%3 == 2)
      v /= 3.21;
  }
  return v;
}

int main(int argc, char *argv[])
{
  LOG(DEBUG) << "debug" << f(1);
  LOG(INFO) << "info" << f(2);
  VLOG(1) << "vlog1 " << f(3);
  LOG(ERROR) << "logerror" << f(4);
  
  cout << "double vector size: " << Vc::double_v::size() << ", float vector size: " << Vc::float_v::size() << endl;

  // test Vc
  int nFiberPoints = 100;
  int nVcVectors = (nFiberPoints + Vc::double_v::Size - 1) / Vc::double_v::Size;

  std::vector<FiberPointData> values(nVcVectors);
  
  for (int i = 0; i < nVcVectors; i++)
  {
    values[i].values[0] = 8.0;
  }

  //Vc::Memory<FiberPointData> data(nFiberPoints);
  
  // test SEMT 
  /*
  // Define an instance of variable with index 0 named x.
  DVAR(x, 0);

  //auto f = pow(x, INT(1)/INT(2));
  //auto f = ln(pow(x, INT(2)));
  auto f = pow(sqrt(x), INT(2));
  cout << "f = " << f << endl;
  
  auto df1 = deriv_t(f, x);
  cout << "f' = " << df1 << endl;

  auto g = x - (x + x);
  cout << "g = " << g << endl;
  
  auto dg1 = deriv_t(g, x);
  cout << "g' = " << dg1 << endl;

  */
  /*
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


  test->debug();*/
  
  return EXIT_SUCCESS;
}
