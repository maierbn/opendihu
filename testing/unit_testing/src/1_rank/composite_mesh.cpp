#include <Python.h>  // this has to be the first included header

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <cassert>
#include <cmath>

#include "gtest/gtest.h"
#include "opendihu.h"
#include "arg.h"
#include "../utility.h"
#include "mesh/face_t.h"
#include "function_space/function_space.h"

TEST(CompositeTest, TwoMeshesWorks0)
{
  // run serial problem
  std::string pythonConfig = R"(

config = {
  "Meshes": {
    "submesh0": {
      "nElements": [3, 2],
      "inputMeshIsGlobal": True,
      "physicalExtent": [2.0, 1.0],
    },
    "submesh1": {
      "nElements": [2, 2],
      "inputMeshIsGlobal": True,
      "physicalExtent": [2.0, 1.0],
      "physicalOffset": [2.0, 1.0],
    },
  },
  "FiniteElementMethod": {
    "inputMeshIsGlobal": True,

    "meshName": ["submesh0", "submesh1"],
  },
}
)";

  DihuContext settings(argc, argv, pythonConfig);

  typedef SpatialDiscretization::FiniteElementMethod<
    Mesh::CompositeOfDimension<2>,
    BasisFunction::LagrangeOfOrder<2>,
    Quadrature::Gauss<1>,
    Equation::None
  > ProblemType;
  ProblemType problem(settings);

  problem.initialize();

  std::string stringRespresentation = problem.data().functionSpace()->meshPartition()->getString();
  std::cout << "stringRespresentation: \n" << stringRespresentation;

  std::string reference = "CompositeMesh, nSubMeshes: 2, removedSharedNodes: [, [1,0] -> [0,34] , ], nElementsLocal: 10, nElementsGlobal: 10, elementNoGlobalBegin: 0, nNodesSharedLocal: 1, nGhostNodesSharedLocal: 0, nRemovedNodesNonGhost: [0, 1, ], nNonDuplicateNodesWithoutGhosts: [35, 24, ], nNodesLocalWithoutGhosts: 59, nNodesLocalWithGhosts: 59, nNodesGlobal: 59, nonDuplicateNodeNoGlobalBegin: 0, meshAndNodeNoLocalToNodeNoNonDuplicateGlobal: [[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,],[-1,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,],], meshAndNodeNoLocalToNodeNoNonDuplicateLocal: [[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,],[-1,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,],], isDuplicate: [[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,],[1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,],], nodeNoNonDuplicateLocalToMeshAndDuplicateLocal: [<0,0>,<0,1>,<0,2>,<0,3>,<0,4>,<0,5>,<0,6>,<0,7>,<0,8>,<0,9>,<0,10>,<0,11>,<0,12>,<0,13>,<0,14>,<0,15>,<0,16>,<0,17>,<0,18>,<0,19>,<0,20>,<0,21>,<0,22>,<0,23>,<0,24>,<0,25>,<0,26>,<0,27>,<0,28>,<0,29>,<0,30>,<0,31>,<0,32>,<0,33>,<0,34>,<1,1>,<1,2>,<1,3>,<1,4>,<1,5>,<1,6>,<1,7>,<1,8>,<1,9>,<1,10>,<1,11>,<1,12>,<1,13>,<1,14>,<1,15>,<1,16>,<1,17>,<1,18>,<1,19>,<1,20>,<1,21>,<1,22>,<1,23>,<1,24>,], nonDuplicateGhostNodeNosGlobal: [], onlyNodalDofLocalNos: [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,], ghostDofNosGlobalPetsc: [], nElementsLocal(): 10, nElementsGlobal(): 10, nDofsLocalWithGhosts(): 59, nDofsLocalWithoutGhosts(): 59, nDofsGlobal(): 59, nNodesLocalWithGhosts(): 59, nNodesLocalWithoutGhosts(): 59, nNodesGlobal(): 59, beginNodeGlobalPetsc(): 0, dofNosLocal(true): [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,], dofNosLocal(false): [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,], ghostDofNosGlobalPetsc(): [], getElementNoGlobalNatural: 0->0 1->1 2->2 3->3 4->4 5->5 6->6 7->7 8->8 9->9 , getNodeNoGlobalPetsc: 0->0 1->1 2->2 3->3 4->4 5->5 6->6 7->7 8->8 9->9 10->10 11->11 12->12 13->13 14->14 15->15 16->16 17->17 18->18 19->19 20->20 21->21 22->22 23->23 24->24 25->25 26->26 27->27 28->28 29->29 30->30 31->31 32->32 33->33 34->34 35->35 36->36 37->37 38->38 39->39 40->40 41->41 42->42 43->43 44->44 45->45 46->46 47->47 48->48 49->49 50->50 51->51 52->52 53->53 54->54 55->55 56->56 57->57 58->58 , getNodeNoGlobalNatural: 0,0->0 0,1->1 0,2->2 0,3->7 0,4->8 0,5->9 0,6->14 0,7->15 0,8->16 1,0->2 1,1->3 1,2->4 1,3->9 1,4->10 1,5->11 1,6->16 1,7->17 1,8->18 2,0->4 2,1->5 2,2->6 2,3->11 2,4->12 2,5->13 2,6->18 2,7->19 2,8->20 3,0->14 3,1->15 3,2->16 3,3->21 3,4->22 3,5->23 3,6->28 3,7->29 3,8->30 4,0->16 4,1->17 4,2->18 4,3->23 4,4->24 4,5->25 4,6->30 4,7->31 4,8->32 5,0->18 5,1->19 5,2->20 5,3->25 5,4->26 5,5->27 5,6->32 5,7->33 5,8->34 6,0->35 6,1->36 6,2->37 6,3->40 6,4->41 6,5->42 6,6->45 6,7->46 6,8->47 7,0->37 7,1->38 7,2->39 7,3->42 7,4->43 7,5->44 7,6->47 7,7->48 7,8->49 8,0->45 8,1->46 8,2->47 8,3->50 8,4->51 8,5->52 8,6->55 8,7->56 8,8->57 9,0->47 9,1->48 9,2->49 9,3->52 9,4->53 9,5->54 9,6->57 9,7->58 9,8->59 , getDofNoGlobalPetsc: [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,], getDofNoGlobalPetsc: 0->0 1->1 2->2 3->3 4->4 5->5 6->6 7->7 8->8 9->9 10->10 11->11 12->12 13->13 14->14 15->15 16->16 17->17 18->18 19->19 20->20 21->21 22->22 23->23 24->24 25->25 26->26 27->27 28->28 29->29 30->30 31->31 32->32 33->33 34->34 35->35 36->36 37->37 38->38 39->39 40->40 41->41 42->42 43->43 44->44 45->45 46->46 47->47 48->48 49->49 50->50 51->51 52->52 53->53 54->54 55->55 56->56 57->57 58->58 , getElementNoLocal: 0->0->1 1->1->1 2->2->1 3->3->1 4->4->1 5->5->1 6->6->1 7->7->1 8->8->1 9->9->1 , getNodeNoLocal: 0->0->1 1->1->1 2->2->1 3->3->1 4->4->1 5->5->1 6->6->1 7->7->1 8->8->1 9->9->1 10->10->1 11->11->1 12->12->1 13->13->1 14->14->1 15->15->1 16->16->1 17->17->1 18->18->1 19->19->1 20->20->1 21->21->1 22->22->1 23->23->1 24->24->1 25->25->1 26->26->1 27->27->1 28->28->1 29->29->1 30->30->1 31->31->1 32->32->1 33->33->1 34->34->1 35->35->1 36->36->1 37->37->1 38->38->1 39->39->1 40->40->1 41->41->1 42->42->1 43->43->1 44->44->1 45->45->1 46->46->1 47->47->1 48->48->1 49->49->1 50->50->1 51->51->1 52->52->1 53->53->1 54->54->1 55->55->1 56->56->1 57->57->1 58->58->1 , getDofNoLocal: 0->0->1 1->1->1 2->2->1 3->3->1 4->4->1 5->5->1 6->6->1 7->7->1 8->8->1 9->9->1 10->10->1 11->11->1 12->12->1 13->13->1 14->14->1 15->15->1 16->16->1 17->17->1 18->18->1 19->19->1 20->20->1 21->21->1 22->22->1 23->23->1 24->24->1 25->25->1 26->26->1 27->27->1 28->28->1 29->29->1 30->30->1 31->31->1 32->32->1 33->33->1 34->34->1 35->35->1 36->36->1 37->37->1 38->38->1 39->39->1 40->40->1 41->41->1 42->42->1 43->43->1 44->44->1 45->45->1 46->46->1 47->47->1 48->48->1 49->49->1 50->50->1 51->51->1 52->52->1 53->53->1 54->54->1 55->55->1 56->56->1 57->57->1 58->58->1 , extractLocalNodesWithoutGhosts: [0,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,25,26,27,28,29,30,31,32,33,34,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,], extractLocalDofsWithoutGhosts: [0,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,25,26,27,28,29,30,31,32,33,34,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,], getSubMeshNoAndElementNoLocal: 0->0,0 1->0,1 2->0,2 3->0,3 4->0,4 5->0,5 6->1,0 7->1,1 8->1,2 9->1,3 , getSubMeshNoAndNodeNoLocal: 0->0,0 1->0,1 2->0,2 3->0,3 4->0,4 5->0,5 6->0,6 7->0,7 8->0,8 9->0,9 10->0,10 11->0,11 12->0,12 13->0,13 14->0,14 15->0,15 16->0,16 17->0,17 18->0,18 19->0,19 20->0,20 21->0,21 22->0,22 23->0,23 24->0,24 25->0,25 26->0,26 27->0,27 28->0,28 29->0,29 30->0,30 31->0,31 32->0,32 33->0,33 34->0,34 35->1,1 36->1,2 37->1,3 38->1,4 39->1,5 40->1,6 41->1,7 42->1,8 43->1,9 44->1,10 45->1,11 46->1,12 47->1,13 48->1,14 49->1,15 50->1,16 51->1,17 52->1,18 53->1,19 54->1,20 55->1,21 56->1,22 57->1,23 58->1,24 , getSubMeshesWithNodes: 0->[<0,0> ] 1->[<0,1> ] 2->[<0,2> ] 3->[<0,3> ] 4->[<0,4> ] 5->[<0,5> ] 6->[<0,6> ] 7->[<0,7> ] 8->[<0,8> ] 9->[<0,9> ] 10->[<0,10> ] 11->[<0,11> ] 12->[<0,12> ] 13->[<0,13> ] 14->[<0,14> ] 15->[<0,15> ] 16->[<0,16> ] 17->[<0,17> ] 18->[<0,18> ] 19->[<0,19> ] 20->[<0,20> ] 21->[<0,21> ] 22->[<0,22> ] 23->[<0,23> ] 24->[<0,24> ] 25->[<0,25> ] 26->[<0,26> ] 27->[<0,27> ] 28->[<0,28> ] 29->[<0,29> ] 30->[<0,30> ] 31->[<0,31> ] 32->[<0,32> ] 33->[<0,33> ] 34->[<0,34> <1,0> ] 35->[<1,1> ] 36->[<1,2> ] 37->[<1,3> ] 38->[<1,4> ] 39->[<1,5> ] 40->[<1,6> ] 41->[<1,7> ] 42->[<1,8> ] 43->[<1,9> ] 44->[<1,10> ] 45->[<1,11> ] 46->[<1,12> ] 47->[<1,13> ] 48->[<1,14> ] 49->[<1,15> ] 50->[<1,16> ] 51->[<1,17> ] 52->[<1,18> ] 53->[<1,19> ] 54->[<1,20> ] 55->[<1,21> ] 56->[<1,22> ] 57->[<1,23> ] 58->[<1,24> ] ";

  ASSERT_EQ(stringRespresentation, reference);
}

TEST(CompositeTest, TwoMeshesWorks1)
{
  // run serial problem
  std::string pythonConfig = R"(

config = {
  "Meshes": {
    "submesh0": {
      "nElements": [3, 2],
      "inputMeshIsGlobal": True,
      "physicalExtent": [2.0, 1.0],
    },
    "submesh1": {
      "nElements": [2, 2],
      "inputMeshIsGlobal": True,
      "physicalExtent": [2.0, 1.0],
      "physicalOffset": [2.0, 0.5],
    },
  },
  "FiniteElementMethod": {
    "inputMeshIsGlobal": True,

    "meshName": ["submesh0", "submesh1"],
  },
}
)";

  DihuContext settings(argc, argv, pythonConfig);

  typedef SpatialDiscretization::FiniteElementMethod<
    Mesh::CompositeOfDimension<2>,
    BasisFunction::LagrangeOfOrder<2>,
    Quadrature::Gauss<1>,
    Equation::None
  > ProblemType;
  ProblemType problem(settings);

  problem.initialize();

  std::string stringRespresentation = problem.data().functionSpace()->meshPartition()->getString();
  std::cout << "stringRespresentation: \n" << stringRespresentation;

  std::string reference = "CompositeMesh, nSubMeshes: 2, removedSharedNodes: [, [1,0] -> [0,20] [1,5] -> [0,27] [1,10] -> [0,34] , ], nElementsLocal: 10, nElementsGlobal: 10, elementNoGlobalBegin: 0, nNodesSharedLocal: 3, nGhostNodesSharedLocal: 0, nRemovedNodesNonGhost: [0, 3, ], nNonDuplicateNodesWithoutGhosts: [35, 22, ], nNodesLocalWithoutGhosts: 57, nNodesLocalWithGhosts: 57, nNodesGlobal: 57, nonDuplicateNodeNoGlobalBegin: 0, meshAndNodeNoLocalToNodeNoNonDuplicateGlobal: [[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,],[-1,35,36,37,38,-1,39,40,41,42,-1,43,44,45,46,47,48,49,50,51,52,53,54,55,56,],], meshAndNodeNoLocalToNodeNoNonDuplicateLocal: [[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,],[-1,35,36,37,38,-1,39,40,41,42,-1,43,44,45,46,47,48,49,50,51,52,53,54,55,56,],], isDuplicate: [[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,],[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,],], nodeNoNonDuplicateLocalToMeshAndDuplicateLocal: [<0,0>,<0,1>,<0,2>,<0,3>,<0,4>,<0,5>,<0,6>,<0,7>,<0,8>,<0,9>,<0,10>,<0,11>,<0,12>,<0,13>,<0,14>,<0,15>,<0,16>,<0,17>,<0,18>,<0,19>,<0,20>,<0,21>,<0,22>,<0,23>,<0,24>,<0,25>,<0,26>,<0,27>,<0,28>,<0,29>,<0,30>,<0,31>,<0,32>,<0,33>,<0,34>,<1,1>,<1,2>,<1,3>,<1,4>,<1,6>,<1,7>,<1,8>,<1,9>,<1,11>,<1,12>,<1,13>,<1,14>,<1,15>,<1,16>,<1,17>,<1,18>,<1,19>,<1,20>,<1,21>,<1,22>,<1,23>,<1,24>,], nonDuplicateGhostNodeNosGlobal: [], onlyNodalDofLocalNos: [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,], ghostDofNosGlobalPetsc: [], nElementsLocal(): 10, nElementsGlobal(): 10, nDofsLocalWithGhosts(): 57, nDofsLocalWithoutGhosts(): 57, nDofsGlobal(): 57, nNodesLocalWithGhosts(): 57, nNodesLocalWithoutGhosts(): 57, nNodesGlobal(): 57, beginNodeGlobalPetsc(): 0, dofNosLocal(true): [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,], dofNosLocal(false): [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,], ghostDofNosGlobalPetsc(): [], getElementNoGlobalNatural: 0->0 1->1 2->2 3->3 4->4 5->5 6->6 7->7 8->8 9->9 , getNodeNoGlobalPetsc: 0->0 1->1 2->2 3->3 4->4 5->5 6->6 7->7 8->8 9->9 10->10 11->11 12->12 13->13 14->14 15->15 16->16 17->17 18->18 19->19 20->20 21->21 22->22 23->23 24->24 25->25 26->26 27->27 28->28 29->29 30->30 31->31 32->32 33->33 34->34 35->35 36->36 37->37 38->38 39->39 40->40 41->41 42->42 43->43 44->44 45->45 46->46 47->47 48->48 49->49 50->50 51->51 52->52 53->53 54->54 55->55 56->56 , getNodeNoGlobalNatural: 0,0->0 0,1->1 0,2->2 0,3->7 0,4->8 0,5->9 0,6->14 0,7->15 0,8->16 1,0->2 1,1->3 1,2->4 1,3->9 1,4->10 1,5->11 1,6->16 1,7->17 1,8->18 2,0->4 2,1->5 2,2->6 2,3->11 2,4->12 2,5->13 2,6->18 2,7->19 2,8->20 3,0->14 3,1->15 3,2->16 3,3->21 3,4->22 3,5->23 3,6->28 3,7->29 3,8->30 4,0->16 4,1->17 4,2->18 4,3->23 4,4->24 4,5->25 4,6->30 4,7->31 4,8->32 5,0->18 5,1->19 5,2->20 5,3->25 5,4->26 5,5->27 5,6->32 5,7->33 5,8->34 6,0->35 6,1->36 6,2->37 6,3->40 6,4->41 6,5->42 6,6->45 6,7->46 6,8->47 7,0->37 7,1->38 7,2->39 7,3->42 7,4->43 7,5->44 7,6->47 7,7->48 7,8->49 8,0->45 8,1->46 8,2->47 8,3->50 8,4->51 8,5->52 8,6->55 8,7->56 8,8->57 9,0->47 9,1->48 9,2->49 9,3->52 9,4->53 9,5->54 9,6->57 9,7->58 9,8->59 , getDofNoGlobalPetsc: [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,], getDofNoGlobalPetsc: 0->0 1->1 2->2 3->3 4->4 5->5 6->6 7->7 8->8 9->9 10->10 11->11 12->12 13->13 14->14 15->15 16->16 17->17 18->18 19->19 20->20 21->21 22->22 23->23 24->24 25->25 26->26 27->27 28->28 29->29 30->30 31->31 32->32 33->33 34->34 35->35 36->36 37->37 38->38 39->39 40->40 41->41 42->42 43->43 44->44 45->45 46->46 47->47 48->48 49->49 50->50 51->51 52->52 53->53 54->54 55->55 56->56 , getElementNoLocal: 0->0->1 1->1->1 2->2->1 3->3->1 4->4->1 5->5->1 6->6->1 7->7->1 8->8->1 9->9->1 , getNodeNoLocal: 0->0->1 1->1->1 2->2->1 3->3->1 4->4->1 5->5->1 6->6->1 7->7->1 8->8->1 9->9->1 10->10->1 11->11->1 12->12->1 13->13->1 14->14->1 15->15->1 16->16->1 17->17->1 18->18->1 19->19->1 20->20->1 21->21->1 22->22->1 23->23->1 24->24->1 25->25->1 26->26->1 27->27->1 28->28->1 29->29->1 30->30->1 31->31->1 32->32->1 33->33->1 34->34->1 35->35->1 36->36->1 37->37->1 38->38->1 39->39->1 40->40->1 41->41->1 42->42->1 43->43->1 44->44->1 45->45->1 46->46->1 47->47->1 48->48->1 49->49->1 50->50->1 51->51->1 52->52->1 53->53->1 54->54->1 55->55->1 56->56->1 , getDofNoLocal: 0->0->1 1->1->1 2->2->1 3->3->1 4->4->1 5->5->1 6->6->1 7->7->1 8->8->1 9->9->1 10->10->1 11->11->1 12->12->1 13->13->1 14->14->1 15->15->1 16->16->1 17->17->1 18->18->1 19->19->1 20->20->1 21->21->1 22->22->1 23->23->1 24->24->1 25->25->1 26->26->1 27->27->1 28->28->1 29->29->1 30->30->1 31->31->1 32->32->1 33->33->1 34->34->1 35->35->1 36->36->1 37->37->1 38->38->1 39->39->1 40->40->1 41->41->1 42->42->1 43->43->1 44->44->1 45->45->1 46->46->1 47->47->1 48->48->1 49->49->1 50->50->1 51->51->1 52->52->1 53->53->1 54->54->1 55->55->1 56->56->1 , extractLocalNodesWithoutGhosts: [0,35,36,37,38,5,39,40,41,42,10,43,44,45,46,47,48,49,50,51,52,53,54,55,56,25,26,27,28,29,30,31,32,33,34,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,], extractLocalDofsWithoutGhosts: [0,35,36,37,38,5,39,40,41,42,10,43,44,45,46,47,48,49,50,51,52,53,54,55,56,25,26,27,28,29,30,31,32,33,34,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,], getSubMeshNoAndElementNoLocal: 0->0,0 1->0,1 2->0,2 3->0,3 4->0,4 5->0,5 6->1,0 7->1,1 8->1,2 9->1,3 , getSubMeshNoAndNodeNoLocal: 0->0,0 1->0,1 2->0,2 3->0,3 4->0,4 5->0,5 6->0,6 7->0,7 8->0,8 9->0,9 10->0,10 11->0,11 12->0,12 13->0,13 14->0,14 15->0,15 16->0,16 17->0,17 18->0,18 19->0,19 20->0,20 21->0,21 22->0,22 23->0,23 24->0,24 25->0,25 26->0,26 27->0,27 28->0,28 29->0,29 30->0,30 31->0,31 32->0,32 33->0,33 34->0,34 35->1,1 36->1,2 37->1,3 38->1,4 39->1,6 40->1,7 41->1,8 42->1,9 43->1,11 44->1,12 45->1,13 46->1,14 47->1,15 48->1,16 49->1,17 50->1,18 51->1,19 52->1,20 53->1,21 54->1,22 55->1,23 56->1,24 , getSubMeshesWithNodes: 0->[<0,0> ] 1->[<0,1> ] 2->[<0,2> ] 3->[<0,3> ] 4->[<0,4> ] 5->[<0,5> ] 6->[<0,6> ] 7->[<0,7> ] 8->[<0,8> ] 9->[<0,9> ] 10->[<0,10> ] 11->[<0,11> ] 12->[<0,12> ] 13->[<0,13> ] 14->[<0,14> ] 15->[<0,15> ] 16->[<0,16> ] 17->[<0,17> ] 18->[<0,18> ] 19->[<0,19> ] 20->[<0,20> <1,0> ] 21->[<0,21> ] 22->[<0,22> ] 23->[<0,23> ] 24->[<0,24> ] 25->[<0,25> ] 26->[<0,26> ] 27->[<0,27> <1,5> ] 28->[<0,28> ] 29->[<0,29> ] 30->[<0,30> ] 31->[<0,31> ] 32->[<0,32> ] 33->[<0,33> ] 34->[<0,34> <1,10> ] 35->[<1,1> ] 36->[<1,2> ] 37->[<1,3> ] 38->[<1,4> ] 39->[<1,6> ] 40->[<1,7> ] 41->[<1,8> ] 42->[<1,9> ] 43->[<1,11> ] 44->[<1,12> ] 45->[<1,13> ] 46->[<1,14> ] 47->[<1,15> ] 48->[<1,16> ] 49->[<1,17> ] 50->[<1,18> ] 51->[<1,19> ] 52->[<1,20> ] 53->[<1,21> ] 54->[<1,22> ] 55->[<1,23> ] 56->[<1,24> ] ";

  ASSERT_EQ(stringRespresentation, reference);
}

TEST(CompositeTest, ThreeMeshesWorks0)
{
  // run serial problem
  std::string pythonConfig = R"(

config = {
  "Meshes": {
    "submesh0": {
      "nElements": [2, 1],
      "inputMeshIsGlobal": True,
      "physicalExtent": [2.0, 1.0],
    },
    "submesh1": {
      "nElements": [1, 1],
      "inputMeshIsGlobal": True,
      "physicalExtent": [1.0, 1.0],
      "physicalOffset": [2.0, 0.5],
    },
    "submesh2": {
      "nElements": [2, 2],
      "inputMeshIsGlobal": True,
      "physicalExtent": [2.0, 2.0],
      "physicalOffset": [2.0, -1.5],
    },
  },
  "FiniteElementMethod": {
    "inputMeshIsGlobal": True,

    "meshName": ["submesh0", "submesh1", "submesh2"],
  },
}
)";

  DihuContext settings(argc, argv, pythonConfig);

  typedef SpatialDiscretization::FiniteElementMethod<
    Mesh::CompositeOfDimension<2>,
    BasisFunction::LagrangeOfOrder<2>,
    Quadrature::Gauss<1>,
    Equation::None
  > ProblemType;
  ProblemType problem(settings);

  problem.initialize();

  std::string stringRespresentation = problem.data().functionSpace()->meshPartition()->getString();
  std::cout << "stringRespresentation: \n" << stringRespresentation;

  std::string reference = "CompositeMesh, nSubMeshes: 3, removedSharedNodes: [, [1,0] -> [0,9] [1,3] -> [0,14] , [2,15] -> [0,4] [2,20] -> [0,9] [2,21] -> [1,1] [2,22] -> [1,2] , ], nElementsLocal: 7, nElementsGlobal: 7, elementNoGlobalBegin: 0, nNodesSharedLocal: 6, nGhostNodesSharedLocal: 0, nRemovedNodesNonGhost: [0, 2, 4, ], nNonDuplicateNodesWithoutGhosts: [15, 7, 21, ], nNodesLocalWithoutGhosts: 43, nNodesLocalWithGhosts: 43, nNodesGlobal: 43, nonDuplicateNodeNoGlobalBegin: 0, meshAndNodeNoLocalToNodeNoNonDuplicateGlobal: [[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,],[-1,15,16,-1,17,18,19,20,21,],[22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,-1,37,38,39,40,-1,-1,-1,41,42,],], meshAndNodeNoLocalToNodeNoNonDuplicateLocal: [[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,],[-1,15,16,-1,17,18,19,20,21,],[22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,-1,37,38,39,40,-1,-1,-1,41,42,],], isDuplicate: [[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,],[1,0,0,1,0,0,0,0,0,],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,1,1,0,0,],], nodeNoNonDuplicateLocalToMeshAndDuplicateLocal: [<0,0>,<0,1>,<0,2>,<0,3>,<0,4>,<0,5>,<0,6>,<0,7>,<0,8>,<0,9>,<0,10>,<0,11>,<0,12>,<0,13>,<0,14>,<1,1>,<1,2>,<1,4>,<1,5>,<1,6>,<1,7>,<1,8>,<2,0>,<2,1>,<2,2>,<2,3>,<2,4>,<2,5>,<2,6>,<2,7>,<2,8>,<2,9>,<2,10>,<2,11>,<2,12>,<2,13>,<2,14>,<2,16>,<2,17>,<2,18>,<2,19>,<2,23>,<2,24>,], nonDuplicateGhostNodeNosGlobal: [], onlyNodalDofLocalNos: [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,], ghostDofNosGlobalPetsc: [], nElementsLocal(): 7, nElementsGlobal(): 7, nDofsLocalWithGhosts(): 43, nDofsLocalWithoutGhosts(): 43, nDofsGlobal(): 43, nNodesLocalWithGhosts(): 43, nNodesLocalWithoutGhosts(): 43, nNodesGlobal(): 43, beginNodeGlobalPetsc(): 0, dofNosLocal(true): [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,], dofNosLocal(false): [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,], ghostDofNosGlobalPetsc(): [], getElementNoGlobalNatural: 0->0 1->1 2->2 3->3 4->4 5->5 6->6 , getNodeNoGlobalPetsc: 0->0 1->1 2->2 3->3 4->4 5->5 6->6 7->7 8->8 9->9 10->10 11->11 12->12 13->13 14->14 15->15 16->16 17->17 18->18 19->19 20->20 21->21 22->22 23->23 24->24 25->25 26->26 27->27 28->28 29->29 30->30 31->31 32->32 33->33 34->34 35->35 36->36 37->37 38->38 39->39 40->40 41->41 42->42 , getNodeNoGlobalNatural: 0,0->0 0,1->1 0,2->2 0,3->5 0,4->6 0,5->7 0,6->10 0,7->11 0,8->12 1,0->2 1,1->3 1,2->4 1,3->7 1,4->8 1,5->9 1,6->12 1,7->13 1,8->14 2,0->15 2,1->16 2,2->17 2,3->18 2,4->19 2,5->20 2,6->21 2,7->22 2,8->23 3,0->24 3,1->25 3,2->26 3,3->29 3,4->30 3,5->31 3,6->34 3,7->35 3,8->36 4,0->26 4,1->27 4,2->28 4,3->31 4,4->32 4,5->33 4,6->36 4,7->37 4,8->38 5,0->34 5,1->35 5,2->36 5,3->39 5,4->40 5,5->41 5,6->44 5,7->45 5,8->46 6,0->36 6,1->37 6,2->38 6,3->41 6,4->42 6,5->43 6,6->46 6,7->47 6,8->48 , getDofNoGlobalPetsc: [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,], getDofNoGlobalPetsc: 0->0 1->1 2->2 3->3 4->4 5->5 6->6 7->7 8->8 9->9 10->10 11->11 12->12 13->13 14->14 15->15 16->16 17->17 18->18 19->19 20->20 21->21 22->22 23->23 24->24 25->25 26->26 27->27 28->28 29->29 30->30 31->31 32->32 33->33 34->34 35->35 36->36 37->37 38->38 39->39 40->40 41->41 42->42 , getElementNoLocal: 0->0->1 1->1->1 2->2->1 3->3->1 4->4->1 5->5->1 6->6->1 , getNodeNoLocal: 0->0->1 1->1->1 2->2->1 3->3->1 4->4->1 5->5->1 6->6->1 7->7->1 8->8->1 9->9->1 10->10->1 11->11->1 12->12->1 13->13->1 14->14->1 15->15->1 16->16->1 17->17->1 18->18->1 19->19->1 20->20->1 21->21->1 22->22->1 23->23->1 24->24->1 25->25->1 26->26->1 27->27->1 28->28->1 29->29->1 30->30->1 31->31->1 32->32->1 33->33->1 34->34->1 35->35->1 36->36->1 37->37->1 38->38->1 39->39->1 40->40->1 41->41->1 42->42->1 , getDofNoLocal: 0->0->1 1->1->1 2->2->1 3->3->1 4->4->1 5->5->1 6->6->1 7->7->1 8->8->1 9->9->1 10->10->1 11->11->1 12->12->1 13->13->1 14->14->1 15->15->1 16->16->1 17->17->1 18->18->1 19->19->1 20->20->1 21->21->1 22->22->1 23->23->1 24->24->1 25->25->1 26->26->1 27->27->1 28->28->1 29->29->1 30->30->1 31->31->1 32->32->1 33->33->1 34->34->1 35->35->1 36->36->1 37->37->1 38->38->1 39->39->1 40->40->1 41->41->1 42->42->1 , extractLocalNodesWithoutGhosts: [22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,0,37,38,39,40,0,0,0,41,42,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,], extractLocalDofsWithoutGhosts: [22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,0,37,38,39,40,0,0,0,41,42,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,], getSubMeshNoAndElementNoLocal: 0->0,0 1->0,1 2->1,0 3->2,0 4->2,1 5->2,2 6->2,3 , getSubMeshNoAndNodeNoLocal: 0->0,0 1->0,1 2->0,2 3->0,3 4->0,4 5->0,5 6->0,6 7->0,7 8->0,8 9->0,9 10->0,10 11->0,11 12->0,12 13->0,13 14->0,14 15->1,1 16->1,2 17->1,4 18->1,5 19->1,6 20->1,7 21->1,8 22->2,0 23->2,1 24->2,2 25->2,3 26->2,4 27->2,5 28->2,6 29->2,7 30->2,8 31->2,9 32->2,10 33->2,11 34->2,12 35->2,13 36->2,14 37->2,16 38->2,17 39->2,18 40->2,19 41->2,23 42->2,24 , getSubMeshesWithNodes: 0->[<0,0> ] 1->[<0,1> ] 2->[<0,2> ] 3->[<0,3> ] 4->[<0,4> <2,15> ] 5->[<0,5> ] 6->[<0,6> ] 7->[<0,7> ] 8->[<0,8> ] 9->[<0,9> <1,0> <2,20> ] 10->[<0,10> ] 11->[<0,11> ] 12->[<0,12> ] 13->[<0,13> ] 14->[<0,14> <1,3> ] 15->[<1,1> <2,21> ] 16->[<1,2> <2,22> ] 17->[<1,4> ] 18->[<1,5> ] 19->[<1,6> ] 20->[<1,7> ] 21->[<1,8> ] 22->[<2,0> ] 23->[<2,1> ] 24->[<2,2> ] 25->[<2,3> ] 26->[<2,4> ] 27->[<2,5> ] 28->[<2,6> ] 29->[<2,7> ] 30->[<2,8> ] 31->[<2,9> ] 32->[<2,10> ] 33->[<2,11> ] 34->[<2,12> ] 35->[<2,13> ] 36->[<2,14> ] 37->[<2,16> ] 38->[<2,17> ] 39->[<2,18> ] 40->[<2,19> ] 41->[<2,23> ] 42->[<2,24> ] ";

  ASSERT_EQ(stringRespresentation, reference);
}

TEST(CompositeTest, ThreeMeshesWorks1)
{
  // run serial problem
  std::string pythonConfig = R"(

config = {
  "Meshes": {
    "submesh0": {
      "nElements": [2, 1],
      "inputMeshIsGlobal": True,
      "physicalExtent": [2.0, 1.0],
    },
    "submesh1": {
      "nElements": [2, 1],
      "inputMeshIsGlobal": True,
      "physicalExtent": [2.0, 1.0],
      "physicalOffset": [-0.5, -1.0],
    },
    "submesh2": {
      "nElements": [2, 1],
      "inputMeshIsGlobal": True,
      "physicalExtent": [2.0, 1.0],
      "physicalOffset": [0.0, -2.0],
    },
  },
  "FiniteElementMethod": {
    "inputMeshIsGlobal": True,

    "meshName": ["submesh0", "submesh1", "submesh2"],
  },
}
)";

  DihuContext settings(argc, argv, pythonConfig);

  typedef SpatialDiscretization::FiniteElementMethod<
    Mesh::CompositeOfDimension<2>,
    BasisFunction::LagrangeOfOrder<2>,
    Quadrature::Gauss<1>,
    Equation::None
  > ProblemType;
  ProblemType problem(settings);

  problem.initialize();

  std::string stringRespresentation = problem.data().functionSpace()->meshPartition()->getString();
  std::cout << "stringRespresentation: \n" << stringRespresentation;

  std::string reference = "CompositeMesh, nSubMeshes: 3, removedSharedNodes: [, [1,11] -> [0,0] [1,12] -> [0,1] [1,13] -> [0,2] [1,14] -> [0,3] , [2,10] -> [1,1] [2,11] -> [1,2] [2,12] -> [1,3] [2,13] -> [1,4] , ], nElementsLocal: 6, nElementsGlobal: 6, elementNoGlobalBegin: 0, nNodesSharedLocal: 8, nGhostNodesSharedLocal: 0, nRemovedNodesNonGhost: [0, 4, 4, ], nNonDuplicateNodesWithoutGhosts: [15, 11, 11, ], nNodesLocalWithoutGhosts: 37, nNodesLocalWithGhosts: 37, nNodesGlobal: 37, nonDuplicateNodeNoGlobalBegin: 0, meshAndNodeNoLocalToNodeNoNonDuplicateGlobal: [[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,],[15,16,17,18,19,20,21,22,23,24,25,-1,-1,-1,-1,],[26,27,28,29,30,31,32,33,34,35,-1,-1,-1,-1,36,],], meshAndNodeNoLocalToNodeNoNonDuplicateLocal: [[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,],[15,16,17,18,19,20,21,22,23,24,25,-1,-1,-1,-1,],[26,27,28,29,30,31,32,33,34,35,-1,-1,-1,-1,36,],], isDuplicate: [[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,],[0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,],[0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,],], nodeNoNonDuplicateLocalToMeshAndDuplicateLocal: [<0,0>,<0,1>,<0,2>,<0,3>,<0,4>,<0,5>,<0,6>,<0,7>,<0,8>,<0,9>,<0,10>,<0,11>,<0,12>,<0,13>,<0,14>,<1,0>,<1,1>,<1,2>,<1,3>,<1,4>,<1,5>,<1,6>,<1,7>,<1,8>,<1,9>,<1,10>,<2,0>,<2,1>,<2,2>,<2,3>,<2,4>,<2,5>,<2,6>,<2,7>,<2,8>,<2,9>,<2,14>,], nonDuplicateGhostNodeNosGlobal: [], onlyNodalDofLocalNos: [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,], ghostDofNosGlobalPetsc: [], nElementsLocal(): 6, nElementsGlobal(): 6, nDofsLocalWithGhosts(): 37, nDofsLocalWithoutGhosts(): 37, nDofsGlobal(): 37, nNodesLocalWithGhosts(): 37, nNodesLocalWithoutGhosts(): 37, nNodesGlobal(): 37, beginNodeGlobalPetsc(): 0, dofNosLocal(true): [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,], dofNosLocal(false): [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,], ghostDofNosGlobalPetsc(): [], getElementNoGlobalNatural: 0->0 1->1 2->2 3->3 4->4 5->5 , getNodeNoGlobalPetsc: 0->0 1->1 2->2 3->3 4->4 5->5 6->6 7->7 8->8 9->9 10->10 11->11 12->12 13->13 14->14 15->15 16->16 17->17 18->18 19->19 20->20 21->21 22->22 23->23 24->24 25->25 26->26 27->27 28->28 29->29 30->30 31->31 32->32 33->33 34->34 35->35 36->36 , getNodeNoGlobalNatural: 0,0->0 0,1->1 0,2->2 0,3->5 0,4->6 0,5->7 0,6->10 0,7->11 0,8->12 1,0->2 1,1->3 1,2->4 1,3->7 1,4->8 1,5->9 1,6->12 1,7->13 1,8->14 2,0->15 2,1->16 2,2->17 2,3->20 2,4->21 2,5->22 2,6->25 2,7->26 2,8->27 3,0->17 3,1->18 3,2->19 3,3->22 3,4->23 3,5->24 3,6->27 3,7->28 3,8->29 4,0->30 4,1->31 4,2->32 4,3->35 4,4->36 4,5->37 4,6->40 4,7->41 4,8->42 5,0->32 5,1->33 5,2->34 5,3->37 5,4->38 5,5->39 5,6->42 5,7->43 5,8->44 , getDofNoGlobalPetsc: [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,], getDofNoGlobalPetsc: 0->0 1->1 2->2 3->3 4->4 5->5 6->6 7->7 8->8 9->9 10->10 11->11 12->12 13->13 14->14 15->15 16->16 17->17 18->18 19->19 20->20 21->21 22->22 23->23 24->24 25->25 26->26 27->27 28->28 29->29 30->30 31->31 32->32 33->33 34->34 35->35 36->36 , getElementNoLocal: 0->0->1 1->1->1 2->2->1 3->3->1 4->4->1 5->5->1 , getNodeNoLocal: 0->0->1 1->1->1 2->2->1 3->3->1 4->4->1 5->5->1 6->6->1 7->7->1 8->8->1 9->9->1 10->10->1 11->11->1 12->12->1 13->13->1 14->14->1 15->15->1 16->16->1 17->17->1 18->18->1 19->19->1 20->20->1 21->21->1 22->22->1 23->23->1 24->24->1 25->25->1 26->26->1 27->27->1 28->28->1 29->29->1 30->30->1 31->31->1 32->32->1 33->33->1 34->34->1 35->35->1 36->36->1 , getDofNoLocal: 0->0->1 1->1->1 2->2->1 3->3->1 4->4->1 5->5->1 6->6->1 7->7->1 8->8->1 9->9->1 10->10->1 11->11->1 12->12->1 13->13->1 14->14->1 15->15->1 16->16->1 17->17->1 18->18->1 19->19->1 20->20->1 21->21->1 22->22->1 23->23->1 24->24->1 25->25->1 26->26->1 27->27->1 28->28->1 29->29->1 30->30->1 31->31->1 32->32->1 33->33->1 34->34->1 35->35->1 36->36->1 , extractLocalNodesWithoutGhosts: [26,27,28,29,30,31,32,33,34,35,25,11,12,13,36,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,], extractLocalDofsWithoutGhosts: [26,27,28,29,30,31,32,33,34,35,25,11,12,13,36,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,], getSubMeshNoAndElementNoLocal: 0->0,0 1->0,1 2->1,0 3->1,1 4->2,0 5->2,1 , getSubMeshNoAndNodeNoLocal: 0->0,0 1->0,1 2->0,2 3->0,3 4->0,4 5->0,5 6->0,6 7->0,7 8->0,8 9->0,9 10->0,10 11->0,11 12->0,12 13->0,13 14->0,14 15->1,0 16->1,1 17->1,2 18->1,3 19->1,4 20->1,5 21->1,6 22->1,7 23->1,8 24->1,9 25->1,10 26->2,0 27->2,1 28->2,2 29->2,3 30->2,4 31->2,5 32->2,6 33->2,7 34->2,8 35->2,9 36->2,14 , getSubMeshesWithNodes: 0->[<0,0> <1,11> ] 1->[<0,1> <1,12> ] 2->[<0,2> <1,13> ] 3->[<0,3> <1,14> ] 4->[<0,4> ] 5->[<0,5> ] 6->[<0,6> ] 7->[<0,7> ] 8->[<0,8> ] 9->[<0,9> ] 10->[<0,10> ] 11->[<0,11> ] 12->[<0,12> ] 13->[<0,13> ] 14->[<0,14> ] 15->[<1,0> ] 16->[<1,1> <2,10> ] 17->[<1,2> <2,11> ] 18->[<1,3> <2,12> ] 19->[<1,4> <2,13> ] 20->[<1,5> ] 21->[<1,6> ] 22->[<1,7> ] 23->[<1,8> ] 24->[<1,9> ] 25->[<1,10> ] 26->[<2,0> ] 27->[<2,1> ] 28->[<2,2> ] 29->[<2,3> ] 30->[<2,4> ] 31->[<2,5> ] 32->[<2,6> ] 33->[<2,7> ] 34->[<2,8> ] 35->[<2,9> ] 36->[<2,14> ] ";

  ASSERT_EQ(stringRespresentation, reference);
}

