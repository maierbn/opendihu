#include <Python.h>  // this has to be the first included header

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <cassert>

#include "gtest/gtest.h"
#include "opendihu.h"
#include "../utility.h"
#include "arg.h"
#include "stiffness_matrix_tester.h"
#include "node_positions_tester.h"

TEST(SolidMechanicsTest, Test3DLinearElasticity)
{
  std::string pythonConfig = R"(
import numpy as np
import sys, os

nx = 5
ny = 4
nz = 3

# boundary conditions (for linear elements)
dirichlet_bc = {}

# left plane
for k in range(0,nz+1):
  for j in range(0,ny+1):
    dirichlet_bc[k*(nx+1)*(ny+1) + j*(nx+1)] = [0.0,None,None]

# front plane
for k in range(0,nz+1):
  for i in range(0,nx+1):
    dirichlet_bc[k*(nx+1)*(ny+1) + i] = [None,0.0,None]

# bottom plane
for j in range(0,ny+1):
  for i in range(0,nx+1):
    dirichlet_bc[j*(nx+1) + i] = [None,None,0.0]

# vertical edge
for k in range(0,nz+1):
  dirichlet_bc[k*(nx+1)*(ny+1)] = [0.0,0.0,None]

# horizontal edge
for i in range(0,nx+1):
  dirichlet_bc[i] = [None,0.0,0.0]

# horizontal edge
for j in range(0,ny+1):
  dirichlet_bc[j*(nx+1)] = [0.0,None,0.0]

# corner
dirichlet_bc[0] = [0.0,0.0,0.0]

neumann_bc = [{"element": k*nx*ny + j*nx + nx-1, "constantVector": [+0.1,+0.6,3.0], "face": "0+"} for k in range(nz) for j in range(ny)]

#dirichlet_bc = {}
#neumann_bc = []

config = {
  "FiniteElementMethod" : {
    "nElements": [nx, ny, nz],
    "inputMeshIsGlobal": True,
    "physicalExtent": [nx, ny, nz],
    "outputInterval": 1.0,
    "prefactor": 1,
    "dirichletBoundaryConditions": dirichlet_bc,
    "neumannBoundaryConditions": neumann_bc,
    "relativeTolerance": 1e-15,
    "solverType": "gmres",
    "preconditionerType": "none",
    "maxIterations": 1e4,
    "bulkModulus": 1.5,
    "shearModulus": 2.0,
    "OutputWriter" : [
      {"format": "Paraview", "outputInterval": 1, "filename": "out/out", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True},
      {"format": "PythonFile", "filename": "out_solid_mechanics3d", "outputInterval": 1, "binary":False, "onlyNodalValues":True},
    ]
  },
}

)";
  DihuContext settings(argc, argv, pythonConfig);
  
  // linear elasticity

  SpatialDiscretization::FiniteElementMethod<
    Mesh::StructuredDeformableOfDimension<3>,
    BasisFunction::LagrangeOfOrder<1>,
    Quadrature::Gauss<3>,
    Equation::Static::LinearElasticity
  > problem(settings);

  problem.run();

  std::string referenceOutput = R"({"meshType": "StructuredDeformable", "dimension": 3, "nElementsGlobal": [5, 4, 3], "nElementsLocal": [5, 4, 3], "beginNodeGlobalNatural": [0, 0, 0], "hasFullNumberOfNodes": [true, true, true], "basisFunction": "Lagrange", "basisOrder": 1, "onlyNodalValues": true, "nRanks": 1, "ownRankNo": 0, "data": [{"name": "geometry", "components": [{"name": "x", "values": [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0]}, {"name": "y", "values": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0]}, {"name": "z", "values": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0]}]}, {"name": "solution", "components": [{"name": "0", "values": [0.0, 0.26352040749518246, 0.5938272050288149, 1.017804873925747, 1.482850669197994, 1.664979439729417, 0.0, 0.24831229762669804, 0.5632836754389753, 0.9716241162463817, 1.415381702525523, 1.5868576075790388, 0.0, 0.20592057179085824, 0.47722790453478797, 0.842217102003177, 1.247310943099748, 1.4046834916047064, 0.0, 0.1480308600680536, 0.35467703709447257, 0.6493172002305352, 1.0025063235469358, 1.146266290455368, 0.0, 0.10047563185541152, 0.24126787117446533, 0.43600868575199275, 0.6544081139811178, 0.7257449629872726, 0.0, 0.19048064783761173, 0.43287016420989866, 0.7550543582619476, 1.098651068086733, 1.2247239927839169, 0.0, 0.17532191224840207, 0.4024156451248284, 0.7089887886035248, 1.0312969181521205, 1.1466903750106097, 0.0, 0.13303334745349357, 0.31653644614725157, 0.5798174621275844, 0.863473213663834, 0.9646854053346716, 0.0, 0.07499833715083257, 0.19371044352023326, 0.3865286059527679, 0.6184476511159657, 0.7061152339320372, 0.0, 0.026243802341574203, 0.07808242800277648, 0.17020924461920708, 0.26598822708919473, 0.27996432843894936, 0.0, 0.018434934359792373, 0.029228013071634553, 0.04067201094616193, 0.07902132454932233, 0.09469714464276469, 0.0, 0.003359996946862066, -0.001084194150765446, -0.005232266556314193, 0.011861467695672399, 0.016828575950076537, 0.0, -0.03867512002350102, -0.08654796891920716, -0.13395457122672946, -0.1555920057443348, -0.1650337943891098, 0.0, -0.09682137674318045, -0.20959312663704033, -0.32776074242266673, -0.40135983840709905, -0.4246625269141226, 0.0, -0.14863676954045973, -0.33087567777652316, -0.551228782538086, -0.7620926457247235, -0.8602375576596021, 0.0, -0.10654972911778643, -0.3402997366560274, -0.8158640238933212, -1.5096144812931644, -1.8982200779913088, 0.0, -0.12174628904399667, -0.3708695349546478, -0.8621991865367706, -1.5775470114309085, -1.9771874331160693, 0.0, -0.16394357373640422, -0.45673280724017723, -0.9916658364657638, -1.7462035019080118, -2.160728347839527, 0.0, -0.22177399549802213, -0.5792027459247792, -1.1851724511758712, -1.992694831458157, -2.4219599864073493, 0.0, -0.27492485542119477, -0.7031409646390167, -1.4117358969291895, -2.3563782736286445, -2.8607120169420357]}, {"name": "1", "values": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.038250544565852355, -0.030861395171500744, -0.0035146548106351075, 0.059669189019066114, 0.1860212980888547, 0.48304237528244304, -0.0685089210182786, -0.05606989785331466, -0.009060862533347336, 0.1029335298724257, 0.33666337532274254, 0.7533022988430976, -0.084590313874762, -0.07053323416494683, -0.016370855910359743, 0.11711806372638683, 0.4002613451830336, 0.9101084461796263, -0.08700410233717089, -0.07360539487057872, -0.021422464459509717, 0.10836807206315746, 0.3877327420833453, 0.919220764299906, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.037517359972199606, -0.030002325507101153, -0.0023019712944422448, 0.06138665366984586, 0.18825646557828618, 0.4854589111291526, -0.06692575320055213, -0.054154026988329426, -0.006213577902475565, 0.10715803270847919, 0.34236279971167227, 0.7596744823085809, -0.08194403941292384, -0.06720425691123572, -0.011101341436589234, 0.1252351716539537, 0.4121163591532377, 0.9251432918932835, -0.08305375075010545, -0.06852003073775866, -0.013122550590281315, 0.12161328382243779, 0.40782856926792843, 0.9522893148550985, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.03607427169519575, -0.028315627262991666, 5.871531463856549e-05, 0.06468614973851092, 0.19254088952959347, 0.49038937614828987, -0.06351655622864438, -0.05001440749575015, -4.4769994623648855e-05, 0.11620871832237979, 0.3544502015960439, 0.7735038984091658, -0.07606249321645876, -0.05969569734945734, 0.001058534187525962, 0.14410485169790668, 0.4384556025212797, 0.9575475029231814, -0.0741928608513048, -0.056736621731050714, 0.007260434631819323, 0.15471118550785815, 0.45527177834917515, 1.018466402558862, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.03613406500381868, -0.028406010710466917, -0.00014739862892098593, 0.0642000925820607, 0.19135841117681465, 0.48876141660613037, -0.0620479178609888, -0.04825295676365835, 0.0024594294518203923, 0.11961236330702454, 0.3585110543024484, 0.778051391520148, -0.07122589984554646, -0.05351006552867593, 0.011131771230972018, 0.1594859727350388, 0.4588770525131827, 0.9812799463407348, -0.06615680930186207, -0.04559190747065901, 0.02778172891132081, 0.18997366281406008, 0.5013681541815822, 1.0689329995725063]}, {"name": "2", "values": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.14975792002467458, -0.1432772718517055, -0.08741924807225201, 0.13922709836913313, 0.7529352318306652, 2.3267681106085356, -0.1494694675809121, -0.1429200905970529, -0.08685933311066409, 0.14012242362670663, 0.7542153172840705, 2.328186907278392, -0.1491633804403576, -0.14248913221431714, -0.0860521947944054, 0.1415825572086384, 0.7566120518837759, 2.331004577876645, -0.1503076332735408, -0.14375549845105215, -0.08763155951123441, 0.1396297176413717, 0.7548542913324461, 2.3310649991828862, -0.15498384664545528, -0.14960580098312, -0.0968339584537666, 0.12550768019154276, 0.7345137194609447, 2.3078042993594545, -0.2407130837175993, -0.23234217634973892, -0.14804968293847148, 0.21636927767645206, 1.2867289851031718, 3.598183338224744, -0.24006888992230058, -0.23156445451341004, -0.14687874726750696, 0.21818191694618805, 1.2894392603010845, 3.601545850503763, -0.23913385998149486, -0.2303610836700196, -0.14487694330552242, 0.22149061753724542, 1.2945947345024509, 3.608161199729829, -0.24055674155813214, -0.23191625188683046, -0.14672227145338426, 0.21939584229385248, 1.2926227780617452, 3.6078207876481914, -0.2469215443781855, -0.24004925708422964, -0.1599345970493098, 0.1988210465198336, 1.2646605147826089, 3.57733591113094, -0.2630725751936658, -0.2538226479425632, -0.16148237754252023, 0.24813068152188542, 1.4497894454944158, 4.065223362635757, -0.2620339881826784, -0.2525984060504189, -0.15971433716878125, 0.25083490236124756, 1.4539339040691086, 4.072247265907334, -0.2601224881745909, -0.2502758657779874, -0.15617941550874223, 0.256328739752249, 1.4622308942950009, 4.084182226264698, -0.2605828568317835, -0.2507215519829516, -0.15645104658395298, 0.25650536185857764, 1.4628703650060684, 4.085819255209441, -0.26605777074101566, -0.2578823338681234, -0.16849213994792048, 0.23761902578785765, 1.4377463733863256, 4.059066368502148]}]}, {"name": "rightHandSide", "components": [{"name": "0", "values": [0.0, 0.0, 0.0, 0.0, 0.0, -0.025, 0.0, 0.0, 0.0, 0.0, 0.0, -0.05, 0.0, 0.0, 0.0, 0.0, 0.0, -0.05, 0.0, 0.0, 0.0, 0.0, 0.0, -0.05, 0.0, 0.0, 0.0, 0.0, 0.0, -0.025000000000000005, 0.0, 0.0, 0.0, 0.0, 0.0, -0.05, 0.0, 0.0, 0.0, 0.0, 0.0, -0.1, 0.0, 0.0, 0.0, 0.0, 0.0, -0.1, 0.0, 0.0, 0.0, 0.0, 0.0, -0.1, 0.0, 0.0, 0.0, 0.0, 0.0, -0.05, 0.0, 0.0, 0.0, 0.0, 0.0, -0.05, 0.0, 0.0, 0.0, 0.0, 0.0, -0.1, 0.0, 0.0, 0.0, 0.0, 0.0, -0.1, 0.0, 0.0, 0.0, 0.0, 0.0, -0.1, 0.0, 0.0, 0.0, 0.0, 0.0, -0.05, 0.0, 0.0, 0.0, 0.0, 0.0, -0.025, 0.0, 0.0, 0.0, 0.0, 0.0, -0.05, 0.0, 0.0, 0.0, 0.0, 0.0, -0.05, 0.0, 0.0, 0.0, 0.0, 0.0, -0.05, 0.0, 0.0, 0.0, 0.0, 0.0, -0.025]}, {"name": "1", "values": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.3, 0.0, 0.0, 0.0, 0.0, 0.0, -0.3, 0.0, 0.0, 0.0, 0.0, 0.0, -0.3, 0.0, 0.0, 0.0, 0.0, 0.0, -0.15, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.6, 0.0, 0.0, 0.0, 0.0, 0.0, -0.6, 0.0, 0.0, 0.0, 0.0, 0.0, -0.6, 0.0, 0.0, 0.0, 0.0, 0.0, -0.3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.6, 0.0, 0.0, 0.0, 0.0, 0.0, -0.6, 0.0, 0.0, 0.0, 0.0, 0.0, -0.6, 0.0, 0.0, 0.0, 0.0, 0.0, -0.3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.3, 0.0, 0.0, 0.0, 0.0, 0.0, -0.3, 0.0, 0.0, 0.0, 0.0, 0.0, -0.3, 0.0, 0.0, 0.0, 0.0, 0.0, -0.15]}, {"name": "2", "values": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.5, 0.0, 0.0, 0.0, 0.0, 0.0, -3.0, 0.0, 0.0, 0.0, 0.0, 0.0, -3.0, 0.0, 0.0, 0.0, 0.0, 0.0, -3.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.5, 0.0, 0.0, 0.0, 0.0, 0.0, -1.5, 0.0, 0.0, 0.0, 0.0, 0.0, -3.0, 0.0, 0.0, 0.0, 0.0, 0.0, -3.0, 0.0, 0.0, 0.0, 0.0, 0.0, -3.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.5, 0.0, 0.0, 0.0, 0.0, 0.0, -0.75, 0.0, 0.0, 0.0, 0.0, 0.0, -1.5, 0.0, 0.0, 0.0, 0.0, 0.0, -1.5, 0.0, 0.0, 0.0, 0.0, 0.0, -1.5, 0.0, 0.0, 0.0, 0.0, 0.0, -0.75]}]}, {"name": "-rhsNeumannBC", "components": [{"name": "0", "values": [0.0, 0.0, 0.0, 0.0, 0.0, 0.025, 0.0, 0.0, 0.0, 0.0, 0.0, 0.05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.025000000000000005, 0.0, 0.0, 0.0, 0.0, 0.0, 0.05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.025, 0.0, 0.0, 0.0, 0.0, 0.0, 0.05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.025]}, {"name": "1", "values": [0.0, 0.0, 0.0, 0.0, 0.0, 0.15, 0.0, 0.0, 0.0, 0.0, 0.0, 0.3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.15, 0.0, 0.0, 0.0, 0.0, 0.0, 0.3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.15, 0.0, 0.0, 0.0, 0.0, 0.0, 0.3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.15]}, {"name": "2", "values": [0.0, 0.0, 0.0, 0.0, 0.0, 0.75, 0.0, 0.0, 0.0, 0.0, 0.0, 1.5, 0.0, 0.0, 0.0, 0.0, 0.0, 1.5, 0.0, 0.0, 0.0, 0.0, 0.0, 1.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.75, 0.0, 0.0, 0.0, 0.0, 0.0, 1.5, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.5, 0.0, 0.0, 0.0, 0.0, 0.0, 1.5, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.75, 0.0, 0.0, 0.0, 0.0, 0.0, 1.5, 0.0, 0.0, 0.0, 0.0, 0.0, 1.5, 0.0, 0.0, 0.0, 0.0, 0.0, 1.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.75]}]}], "timeStepNo": -1, "currentTime": 0.0})";
  assertFileMatchesContent("out_solid_mechanics3d.py", referenceOutput);
}

TEST(SolidMechanicsTest, Test2DLinearElasticity)
{

  std::string pythonConfig = R"(
import numpy as np

nx = 4
ny = 4

# boundary conditions (for linear elements)
dirichlet_bc = {0: [0.0,0.0,0.0]}

for j in range(1,ny+1):
  dirichlet_bc[j*(nx+1)] = [0.0,None,None]

dirichlet_bc[1] = [0.0,0.0]

neumann_bc = [{"element": j*nx+(nx-1), "constantVector": [0.1,+0.2], "face": "0+"} for j in range(ny)]

#dirichlet_bc = {}
#neumann_bc = []

config = {
  "FiniteElementMethod" : {
    "nElements": [nx, ny],
    "inputMeshIsGlobal": True,
    "physicalExtent": [nx, ny],
    "outputInterval": 1.0,
    "prefactor": 1,
    "dirichletBoundaryConditions": dirichlet_bc,
    "neumannBoundaryConditions": neumann_bc,
    "relativeTolerance": 1e-15,
    "solverType": "gmres",
    "preconditionerType": "none",
    "maxIterations": 1e4,
    "bulkModulus": 1.5,
    "shearModulus": 2.0,
    "OutputWriter" : [
      {"format": "Paraview", "outputInterval": 1, "filename": "out/laplace", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True},
      {"format": "PythonFile", "filename": "out_solid_mechanics2d", "outputInterval": 1, "binary":False, "onlyNodalValues":True},
    ]
  },
}

)";
  DihuContext settings(argc, argv, pythonConfig);

  // linear elasticity

  SpatialDiscretization::FiniteElementMethod<
    Mesh::StructuredDeformableOfDimension<2>,
    BasisFunction::LagrangeOfOrder<1>,
    Quadrature::Gauss<3>,
    Equation::Static::LinearElasticity
  > equationDiscretized(settings);

  equationDiscretized.run();

  std::string referenceOutput = R"({"meshType": "StructuredDeformable", "dimension": 2, "nElementsGlobal": [4, 4], "nElementsLocal": [4, 4], "beginNodeGlobalNatural": [0, 0], "hasFullNumberOfNodes": [true, true], "basisFunction": "Lagrange", "basisOrder": 1, "onlyNodalValues": true, "nRanks": 1, "ownRankNo": 0, "data": [{"name": "geometry", "components": [{"name": "x", "values": [0.0, 1.0, 2.0, 3.0, 4.0, 0.0, 1.0, 2.0, 3.0, 4.0, 0.0, 1.0, 2.0, 3.0, 4.0, 0.0, 1.0, 2.0, 3.0, 4.0, 0.0, 1.0, 2.0, 3.0, 4.0]}, {"name": "y", "values": [0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0, 3.0, 3.0, 4.0, 4.0, 4.0, 4.0, 4.0]}, {"name": "z", "values": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]}]}, {"name": "solution", "components": [{"name": "0", "values": [0.0, 0.0, 0.3072540434618576, 0.44415263465279203, 0.5050170741073159, 0.0, 0.06186959477926346, 0.1071255940889371, 0.18256753423417685, 0.22311019785135827, 0.0, -0.000224101928177859, 0.0056388880586356474, 0.018445675166766477, 0.04196913125168837, 0.0, -0.06382118517546173, -0.11402792498690408, -0.13591977435781616, -0.12691916126254674, 0.0, -0.15558463223927396, -0.2842080976926419, -0.3654136578193336, -0.3804529613763943]}, {"name": "1", "values": [0.0, 0.0, 0.3077701183934221, 0.6139347501739595, 0.9644247021103929, 0.07388799914196018, 0.12140709923854197, 0.31409333791842897, 0.6155585019236227, 0.9387034710636611, 0.11929155912845028, 0.1797169949244743, 0.35497799720331635, 0.6218514485017068, 0.9117217348152413, 0.14649375126430095, 0.2053779935868313, 0.37421670025797105, 0.6285913637930131, 0.9286383098631888, 0.15713203746804746, 0.21534421709416524, 0.38104897683477834, 0.6304862210036131, 0.9489081939755176]}]}, {"name": "rightHandSide", "components": [{"name": "0", "values": [0.0, 0.0, 0.0, 0.0, -0.05, 0.0, 0.0, 0.0, 0.0, -0.1, 0.0, 0.0, 0.0, 0.0, -0.1, 0.0, 0.0, 0.0, 0.0, -0.1, 0.0, 0.0, 0.0, 0.0, -0.05]}, {"name": "1", "values": [0.0, 0.0, 0.0, 0.0, -0.1, 0.0, 0.0, 0.0, 0.0, -0.2, 0.0, 0.0, 0.0, 0.0, -0.2, 0.0, 0.0, 0.0, 0.0, -0.2, 0.0, 0.0, 0.0, 0.0, -0.1]}]}, {"name": "-rhsNeumannBC", "components": [{"name": "0", "values": [0.0, 0.0, 0.0, 0.0, 0.05, 0.0, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0, 0.0, 0.05]}, {"name": "1", "values": [0.0, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0, 0.0, 0.2, 0.0, 0.0, 0.0, 0.0, 0.2, 0.0, 0.0, 0.0, 0.0, 0.2, 0.0, 0.0, 0.0, 0.0, 0.1]}]}], "timeStepNo": -1, "currentTime": 0.0})";
  assertFileMatchesContent("out_solid_mechanics2d.py", referenceOutput);
}

TEST(SolidMechanicsTest, TestFEBio1)
{
  // this test corresponds to the example solid_mechanics/mooney_rivlin_febio
  // check that febio and opendihu compute the same displacements for Mooney-Rivlin material

  // --- febio ----
  std::string pythonConfigFebio = R"(

import numpy as np
import sys, os

# parameters
force = 10
material_parameters = [10, 10, 1e6]       # c0, c1, k

# number of elements
nx = 2    # 2
ny = 2    # 2
nz = 5    # 5
physical_extent = [2, 2, 5]


config = {
  "NonlinearElasticitySolverFebio": {
    "durationLogKey": "febio",
    "force": force*(physical_extent[0]*physical_extent[1]),                # factor of force that is applied in axial direction of the muscle
    "materialParameters": material_parameters,   # c0, c1, k for Ψ = c0 * (I1-3) + c1 * (I2-3) + 1/2*k*(log(J))^2
    
    # mesh
    "nElements": [nx, ny, nz],
    "inputMeshIsGlobal": True,
    "physicalExtent": physical_extent,
    "physicalOffset": [0, 0, 0],        # offset/translation where the whole mesh begins
    #"nodePositions": [[i/nx*(1+k/nz), j/ny*(1+k/nz), (k/nz)**1.1*5] for k in range(nz+1) for j in range(ny+1) for i in range(nx+1)],
    
    # output writer
    "OutputWriter" : [
      {"format": "Paraview", "outputInterval": 1, "filename": "out/febio", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
      {"format": "PythonFile", "filename": "out/febio", "outputInterval": 1, "binary":False, "onlyNodalValues":True, "fileNumbering": "incremental"},
    ],
  },
}

)";
  DihuContext settingsFebio(argc, argv, pythonConfigFebio);
  
  // define problem
  TimeSteppingScheme::NonlinearElasticitySolverFebio problemFebio(settingsFebio);
  
  // only run the test if febio is installed 
  if (!problemFebio.isFebioAvailable())
  {
    LOG(INFO) << "febio2 is not installed, abort test.";
    return;
  }

  // run problem
  problemFebio.run();


  // --- opendihu ----
  std::string pythonConfigOpendihu = R"(

# isotropic Mooney Rivlin
import numpy as np
import sys, os

# parameters
force = 10
material_parameters = [10, 10]       # c0, c1

# number of elements
nx = 2    # 2
ny = 2    # 2
nz = 5    # 5

#physical_extent = [1, 1, 5]
physical_extent = [2, 2, 5]

# number of nodes
mx = 2*nx + 1
my = 2*ny + 1
mz = 2*nz + 1

# set the same Dirichlet boundary conditions as the FEBio settings
# boundary conditions (for quadratic elements)
dirichlet_bc = {}

xpos = 0.0
ypos = 0.0
zpos = 0.0

# fix z direction for all
for j in range(0,my):
  for i in range(0,mx):
    k = 0
    dirichlet_bc[k*mx*my + j*(mx) + i] = [None,None,0]

# fix x direction for left row
for j in range(0,my):
  k = 0
  i = 0
  dirichlet_bc[k*mx*my + j*mx + i][0] = 0

# fix y direction for front row
for i in range(0,mx):
  k = 0
  j = 0
  dirichlet_bc[k*mx*my + j*mx + i][1] = 0

# set the Neumann bc's
neumann_bc = [{"element": (nz-1)*nx*ny + j*nx + i, "constantVector": [0,0,force], "face": "2+"} for j in range(ny) for i in range(nx)]

def compare_result(data):
  
  import sys, os
  import py_reader
  import numpy as np

  #print(os.getcwd())
  # load data
  data1 = py_reader.load_data(["out/febio_0000001.py"])
  data2 = py_reader.load_data(["out/opendihu_0000001.py"])
  
  # if files do not exist
  if data1 == [] or data2 == []:
    return

  component_name = "0"
  total_error = 0

  values1 = py_reader.get_values(data1[0], "geometry", component_name)
  values2 = py_reader.get_values(data2[0], "geometry", component_name)
    
  # values2 contains entries for quadratic elements
  # extract the corner values
  n_elements = data2[0]['nElements']
  nx = n_elements[0]
  ny = n_elements[1]
  nz = n_elements[2]
  mx = nx*2 + 1
  my = ny*2 + 1
  mz = nz*2 + 1

  values2_linear = []
  for k in range(nz+1):
    for j in range(ny+1):
      for i in range(nx+1):
        values2_linear.append(values2[2*k * mx * my + 2*j * mx + 2*i])
    
  #print("values1 (febio):    ",list(values1))
  #print("values2 (opendihu): ",values2_linear)
    
  error_rms = np.sqrt(np.mean((values1-values2_linear)**2))
    
  print("rms: {}".format(error_rms))
  with open("rms", "w") as f:
    f.write(str(error_rms))


config = {
  "scenarioName": "3d_box",
  "solverStructureDiagramFile":     "solver_structure.txt",     # output file of a diagram that shows data connection between solvers
  "mappingsBetweenMeshesLogFile":   "out/mappings_between_meshes.txt",  # output file for log of mappings 
  "HyperelasticitySolver": {
    "durationLogKey": "nonlinear",
    
    "materialParameters":         material_parameters,
    "displacementsScalingFactor": 1.0,   # scaling factor for displacements
    "constantBodyForce":          [0.0, 0.0, 0.0],   # body force in whole body region
    "residualNormLogFilename": "log_residual_norm.txt",
    "useAnalyticJacobian": True,
    "useNumericJacobian": False,   # Only works in parallel execution. If both numeric and analytic are enable, it uses the analytic for the preconditioner and the numeric as normal jacobian
      
    "dumpDenseMatlabVariables": False,   # extra output of matlab vectors, x,r, jacobian matrix
    # if useAnalyticJacobian,useNumericJacobian and dumpDenseMatlabVariables are all three true, the analytic and numeric jacobian matrices will get compared to see if there are programming errors for the analytic jacobian
    
    # mesh
    "nElements": [nx, ny, nz],
    "inputMeshIsGlobal": True,
    "physicalExtent": physical_extent,
    "physicalOffset": [0, 0, 0],        # offset/translation where the whole mesh begins
    #"nodePositions": [[(i/(2.*nx))*(1+k/(2.*nz)), (j/(2*nx))*(1+k/(2.*nz)), (k/(2.*nz))**1.1*5] for k in range(mz) for j in range(my) for i in range(mx)],
    
    # nonlinear solver
    "relativeTolerance": 1e-10,         # 1e-10 relative tolerance of the linear solver
    "absoluteTolerance": 1e-10,         # 1e-10 absolute tolerance of the residual of the linear solver    
    "solverType": "preonly",            # type of the linear solver: cg groppcg pipecg pipecgrr cgne nash stcg gltr richardson chebyshev gmres tcqmr fcg pipefcg bcgs ibcgs fbcgs fbcgsr bcgsl cgs tfqmr cr pipecr lsqr preonly qcg bicg fgmres pipefgmres minres symmlq lgmres lcd gcr pipegcr pgmres dgmres tsirm cgls
    "preconditionerType": "lu",         # type of the preconditioner
    "maxIterations": 1e4,               # maximum number of iterations in the linear solver
    "dumpFilename": "",
    "dumpFormat": "matlab",   # default, ascii, matlab
    "snesMaxFunctionEvaluations": 1e8,  # maximum number of function iterations
    "snesMaxIterations": 20,            # maximum number of iterations in the nonlinear solver
    "snesRelativeTolerance": 1e-5,     # relative tolerance of the nonlinear solver
    "snesLineSearchType": "l2",        # type of linesearch, possible values: "bt" "nleqerr" "basic" "l2" "cp" "ncglinear"
    "snesAbsoluteTolerance": 1e-5,     # absolute tolerance of the nonlinear solver
    "snesRebuildJacobianFrequency": 5, # how often the jacobian should be recomputed, -1 indicates NEVER rebuild, 1 means rebuild every time the Jacobian is computed within a single nonlinear solve, 2 means every second time the Jacobian is built etc. -2 means rebuild at next chance but then never again 
    
    #"loadFactors":  [0.1, 0.2, 0.35, 0.47, 0.5, 0.75, 1.0],   # load factors for every timestep
    #"loadFactors": [],                 # no load factors, solve problem directly
    "nNonlinearSolveCalls": 1,         # how often the nonlinear solve should be repeated
    
    # boundary conditions
    "dirichletBoundaryConditions": dirichlet_bc,
    "neumannBoundaryConditions": neumann_bc,
    "divideNeumannBoundaryConditionValuesByTotalArea": True,            # if the given Neumann boundary condition values under "neumannBoundaryConditions" are total forces instead of surface loads and therefore should be scaled by the surface area of all elements where Neumann BC are applied
    #"updateDirichletBoundaryConditionsFunction": update_dirichlet_boundary_conditions,
    "updateDirichletBoundaryConditionsFunction": None,
    "updateDirichletBoundaryConditionsFunctionCallInterval": 1,
    
    "OutputWriter" : [   # output files for displacements function space (quadratic elements)
      {"format": "Paraview", "outputInterval": 1, "filename": "out/opendihu", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
      {"format": "PythonFile", "filename": "out/opendihu", "outputInterval": 1, "binary":False, "onlyNodalValues":True, "fileNumbering": "incremental"},
      {"format": "PythonCallback", "callback": compare_result, "outputInterval": 1},
    ],
    "pressure": None,
    #"LoadIncrements": None,
    #"pressure": {   # output files for pressure function space (linear elements)
    #  "OutputWriter" : [
    #    {"format": "Paraview", "outputInterval": 1, "filename": "out/p", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
    #    {"format": "PythonFile", "filename": "out/p", "outputInterval": 1, "binary":False, "onlyNodalValues":True, "fileNumbering": "incremental"},
    #  ]
    #},
    # output writer for debugging, outputs files after each load increment, the geometry is not changed but u and v are written
    "LoadIncrements": {   
      "OutputWriter" : [
        {"format": "Paraview", "outputInterval": 1, "filename": "out/load_increments", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
      ]
    },
  },
}

)";
  DihuContext settingsOpendihu(argc, argv, pythonConfigOpendihu);
  
  // define problem
  SpatialDiscretization::HyperelasticitySolver<> problemOpendihu(settingsOpendihu);
  
  // run problem
  problemOpendihu.run();


  // parse rms error
  std::ifstream file("rms");
  ASSERT_TRUE(file.is_open());

  std::string line;
  std::getline(file, line);
  double error_rms = atof(line.c_str());

  ASSERT_LE(error_rms, 1e-4);
}
 

TEST(SolidMechanicsTest, TestFEBio2)
{
  // this test corresponds to the example solid_mechanics/mooney_rivlin_febio
  // check that febio and opendihu compute the same displacements for Mooney-Rivlin material

  // --- febio ----
  std::string pythonConfigFebio = R"(
import numpy as np
import sys, os

# parameters
force = 3
material_parameters = [2, 4, 1e6]       # c0, c1, k

# number of elements
nx = 2    # 2
ny = 2    # 2
nz = 5    # 5
physical_extent = [4, 4, 5]


config = {
  "NonlinearElasticitySolverFebio": {
    "durationLogKey": "febio",
    "force": force*(physical_extent[0]*physical_extent[1]),                # factor of force that is applied in axial direction of the muscle
    "materialParameters": material_parameters,   # c0, c1, k for Ψ = c0 * (I1-3) + c1 * (I2-3) + 1/2*k*(log(J))^2
    
    # mesh
    "nElements": [nx, ny, nz],
    "inputMeshIsGlobal": True,
    "physicalExtent": physical_extent,
    "physicalOffset": [0, 0, 0],        # offset/translation where the whole mesh begins
    #"nodePositions": [[i/nx*(1+k/nz), j/ny*(1+k/nz), (k/nz)**1.1*5] for k in range(nz+1) for j in range(ny+1) for i in range(nx+1)],
    
    # output writer
    "OutputWriter" : [
      {"format": "Paraview", "outputInterval": 1, "filename": "out/febio", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
      {"format": "PythonFile", "filename": "out/febio", "outputInterval": 1, "binary":False, "onlyNodalValues":True, "fileNumbering": "incremental"},
    ],
  },
}

)";
  DihuContext settingsFebio(argc, argv, pythonConfigFebio);
  
  // define problem
  TimeSteppingScheme::NonlinearElasticitySolverFebio problemFebio(settingsFebio);
  
  // only run the test if febio is installed 
  if (!problemFebio.isFebioAvailable())
  {
    LOG(INFO) << "febio2 is not installed, abort test.";
    return;
  }

  // run problem
  problemFebio.run();


  // --- opendihu ----
  std::string pythonConfigOpendihu = R"(

# isotropic Mooney Rivlin
import numpy as np
import sys, os

# parameters
force = 3
material_parameters = [2, 4]       # c0, c1

# number of elements
nx = 2    # 2
ny = 2    # 2
nz = 5    # 5

physical_extent = [4, 4, 5]

# number of nodes
mx = 2*nx + 1
my = 2*ny + 1
mz = 2*nz + 1

# set the same Dirichlet boundary conditions as the FEBio settings
# boundary conditions (for quadratic elements)
dirichlet_bc = {}

xpos = 0.0
ypos = 0.0
zpos = 0.0

# fix z direction for all
for j in range(0,my):
  for i in range(0,mx):
    k = 0
    dirichlet_bc[k*mx*my + j*(mx) + i] = [None,None,0]

# fix x direction for left row
for j in range(0,my):
  k = 0
  i = 0
  dirichlet_bc[k*mx*my + j*mx + i][0] = 0

# fix y direction for front row
for i in range(0,mx):
  k = 0
  j = 0
  dirichlet_bc[k*mx*my + j*mx + i][1] = 0

# set the Neumann bc's
neumann_bc = [{"element": (nz-1)*nx*ny + j*nx + i, "constantVector": [0,0,force], "face": "2+"} for j in range(ny) for i in range(nx)]

def compare_result(data):
  
  import sys, os
  import py_reader
  import numpy as np

  #print(os.getcwd())
  # load data
  data1 = py_reader.load_data(["out/febio_0000001.py"])
  data2 = py_reader.load_data(["out/opendihu_0000001.py"])
  
  # if files do not exist
  if data1 == [] or data2 == []:
    return

  component_name = "0"
  total_error = 0

  values1 = py_reader.get_values(data1[0], "geometry", component_name)
  values2 = py_reader.get_values(data2[0], "geometry", component_name)
    
  # values2 contains entries for quadratic elements
  # extract the corner values
  n_elements = data2[0]['nElements']
  nx = n_elements[0]
  ny = n_elements[1]
  nz = n_elements[2]
  mx = nx*2 + 1
  my = ny*2 + 1
  mz = nz*2 + 1

  values2_linear = []
  for k in range(nz+1):
    for j in range(ny+1):
      for i in range(nx+1):
        values2_linear.append(values2[2*k * mx * my + 2*j * mx + 2*i])
    
  #print("values1 (febio):    ",list(values1))
  #print("values2 (opendihu): ",values2_linear)
    
  error_rms = np.sqrt(np.mean((values1-values2_linear)**2))
    
  print("rms: {}".format(error_rms))
  with open("rms", "w") as f:
    f.write(str(error_rms))


config = {
  "scenarioName": "3d_box",
  "solverStructureDiagramFile":     "solver_structure.txt",     # output file of a diagram that shows data connection between solvers
  "mappingsBetweenMeshesLogFile":   "out/mappings_between_meshes.txt",  # output file for log of mappings 
  "HyperelasticitySolver": {
    "durationLogKey": "nonlinear",
    
    "materialParameters":         material_parameters,
    "displacementsScalingFactor": 1.0,   # scaling factor for displacements
    "constantBodyForce":          [0.0, 0.0, 0.0],   # body force in whole body region
    "residualNormLogFilename": "log_residual_norm.txt",
    "useAnalyticJacobian": True,
    "useNumericJacobian": False,   # Only works in parallel execution. If both numeric and analytic are enable, it uses the analytic for the preconditioner and the numeric as normal jacobian
      
    "dumpDenseMatlabVariables": False,   # extra output of matlab vectors, x,r, jacobian matrix
    # if useAnalyticJacobian,useNumericJacobian and dumpDenseMatlabVariables are all three true, the analytic and numeric jacobian matrices will get compared to see if there are programming errors for the analytic jacobian
    
    # mesh
    "nElements": [nx, ny, nz],
    "inputMeshIsGlobal": True,
    "physicalExtent": physical_extent,
    "physicalOffset": [0, 0, 0],        # offset/translation where the whole mesh begins
    #"nodePositions": [[(i/(2.*nx))*(1+k/(2.*nz)), (j/(2*nx))*(1+k/(2.*nz)), (k/(2.*nz))**1.1*5] for k in range(mz) for j in range(my) for i in range(mx)],
    
    # nonlinear solver
    "relativeTolerance": 1e-10,         # 1e-10 relative tolerance of the linear solver
    "absoluteTolerance": 1e-10,         # 1e-10 absolute tolerance of the residual of the linear solver    
    "solverType": "preonly",            # type of the linear solver: cg groppcg pipecg pipecgrr cgne nash stcg gltr richardson chebyshev gmres tcqmr fcg pipefcg bcgs ibcgs fbcgs fbcgsr bcgsl cgs tfqmr cr pipecr lsqr preonly qcg bicg fgmres pipefgmres minres symmlq lgmres lcd gcr pipegcr pgmres dgmres tsirm cgls
    "preconditionerType": "lu",         # type of the preconditioner
    "maxIterations": 1e4,               # maximum number of iterations in the linear solver
    "dumpFilename": "",
    "dumpFormat": "matlab",   # default, ascii, matlab
    "snesMaxFunctionEvaluations": 1e8,  # maximum number of function iterations
    "snesMaxIterations": 20,            # maximum number of iterations in the nonlinear solver
    "snesRelativeTolerance": 1e-5,     # relative tolerance of the nonlinear solver
    "snesLineSearchType": "l2",        # type of linesearch, possible values: "bt" "nleqerr" "basic" "l2" "cp" "ncglinear"
    "snesAbsoluteTolerance": 1e-5,     # absolute tolerance of the nonlinear solver
    "snesRebuildJacobianFrequency": 5, # how often the jacobian should be recomputed, -1 indicates NEVER rebuild, 1 means rebuild every time the Jacobian is computed within a single nonlinear solve, 2 means every second time the Jacobian is built etc. -2 means rebuild at next chance but then never again 
    
    #"loadFactors":  [0.1, 0.2, 0.35, 0.47, 0.5, 0.75, 1.0],   # load factors for every timestep
    #"loadFactors": [],                 # no load factors, solve problem directly
    "nNonlinearSolveCalls": 1,         # how often the nonlinear solve should be repeated
    
    # boundary conditions
    "dirichletBoundaryConditions": dirichlet_bc,
    "neumannBoundaryConditions": neumann_bc,
    "divideNeumannBoundaryConditionValuesByTotalArea": True,            # if the given Neumann boundary condition values under "neumannBoundaryConditions" are total forces instead of surface loads and therefore should be scaled by the surface area of all elements where Neumann BC are applied
    #"updateDirichletBoundaryConditionsFunction": update_dirichlet_boundary_conditions,
    "updateDirichletBoundaryConditionsFunction": None,
    "updateDirichletBoundaryConditionsFunctionCallInterval": 1,
    
    "OutputWriter" : [   # output files for displacements function space (quadratic elements)
      {"format": "Paraview", "outputInterval": 1, "filename": "out/opendihu", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
      {"format": "PythonFile", "filename": "out/opendihu", "outputInterval": 1, "binary":False, "onlyNodalValues":True, "fileNumbering": "incremental"},
      {"format": "PythonCallback", "callback": compare_result, "outputInterval": 1},
    ],
    "pressure": None,
    #"LoadIncrements": None,
    #"pressure": {   # output files for pressure function space (linear elements)
    #  "OutputWriter" : [
    #    {"format": "Paraview", "outputInterval": 1, "filename": "out/p", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
    #    {"format": "PythonFile", "filename": "out/p", "outputInterval": 1, "binary":False, "onlyNodalValues":True, "fileNumbering": "incremental"},
    #  ]
    #},
    # output writer for debugging, outputs files after each load increment, the geometry is not changed but u and v are written
    "LoadIncrements": {   
      "OutputWriter" : [
        {"format": "Paraview", "outputInterval": 1, "filename": "out/load_increments", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
      ]
    },
  },
}

)";
  DihuContext settingsOpendihu(argc, argv, pythonConfigOpendihu);
  
  // define problem
  SpatialDiscretization::HyperelasticitySolver<> problemOpendihu(settingsOpendihu);
  
  // run problem
  problemOpendihu.run();


  // parse rms error
  std::ifstream file("rms");
  ASSERT_TRUE(file.is_open());

  std::string line;
  std::getline(file, line);
  double error_rms = atof(line.c_str());

  ASSERT_LE(error_rms, 1e-4);
}
