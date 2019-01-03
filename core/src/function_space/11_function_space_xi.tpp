#include "function_space/11_function_space_xi.h"

#include <Python.h>  // has to be the first included header
#include <array>
#include "easylogging++.h"
#include "function_space/00_function_space_base_dim.h"
#include <chrono>

namespace FunctionSpace
{

const double POINT_IN_ELEMENT_EPSILON = 1e-4;    // (1e-5 is too small)
const int N_NEWTON_ITERATIONS = 7;    // (4) (5) number of newton iterations to find out if a point is inside an element
const double RESIDUUM_NORM_TOLERANCE = 1e-4;    // usually 1e-2 takes 2-3 iterations
 
// general implementation
template<typename MeshType,typename BasisFunctionType,typename DummyForTraits>
bool FunctionSpaceXi<MeshType,BasisFunctionType,DummyForTraits>::
pointIsInElement(Vec3 point, element_no_t elementNo, std::array<double,MeshType::dim()> &xi)
{
  // timing measurements are disabled, they showed that 'computeApproximateXiForPoint' makes sense and is faster than just initializing the initial guess to 0
#if 0 
  static double durationApproximation = 0.0;
  static double durationNewton = 0.0;
  static int nMeasurements = 0;
 
  auto tStart = std::chrono::steady_clock::now();
#endif

  VLOG(2) << "pointIsInElement(" << point << " element " << elementNo << ")";

  // for 3D mesh and linear Lagrange basis function compute approximate xi by heuristic, else set to 0
  this->computeApproximateXiForPoint(point, elementNo, xi);
   
#if 0 
  auto tEnd = std::chrono::steady_clock::now();
  durationApproximation += std::chrono::duration_cast<std::chrono::duration<double> >(tEnd - tStart).count();
  nMeasurements++;
  
  tStart = std::chrono::steady_clock::now();
#endif

  // define constants
  const int D = MeshType::dim();
  const int nDofsPerElement = FunctionSpaceFunction<MeshType,BasisFunctionType>::nDofsPerElement();
  
  // get geometry field (which are the node positions for Lagrange basis and node positions and derivatives for Hermite)
  std::array<Vec3,nDofsPerElement> geometryValues;
  this->getElementGeometry(elementNo, geometryValues);

  VLOG(2) << "point " << point << ", geometryValues: " << geometryValues;

  Vec3 residuum = point - this->template interpolateValueInElement<3>(geometryValues, xi);
  double residuumNormSquared = MathUtility::normSquared<3>(residuum);
  if (VLOG_IS_ON(2))
  {
    VLOG(2) << " xi0 = " << xi << ", residuum: " << residuum << " (norm: " << sqrt(residuumNormSquared) << ")";
  }
  
  for (int iterationNo = 0; iterationNo < N_NEWTON_ITERATIONS && residuumNormSquared > MathUtility::sqr(RESIDUUM_NORM_TOLERANCE); iterationNo++)
  {
    Tensor2<D> inverseJacobian = this->getInverseJacobian(geometryValues, elementNo, xi);
    xi += inverseJacobian * MathUtility::transformToD<D,3>(residuum);
    
    residuum = point - this->template interpolateValueInElement<3>(geometryValues, xi);
    residuumNormSquared = MathUtility::normSquared<3>(residuum);
    
    if (VLOG_IS_ON(2))
    {
      VLOG(2) << " xi_i = " << xi << ", residuum: " << residuum << " (norm: " << sqrt(residuumNormSquared) << ")";
    }
  }
  
  // Phi(xi) = point
  // Phi(xi) = Phi(xi0) + J*(xi-xi0)  => xi = xi0 + Jinv*(point - Phi(xi0))
  double epsilon = POINT_IN_ELEMENT_EPSILON;
  
  // check if point is inside the element by looking at the value of xi
  bool pointIsInElement = true;
  for (int i = 0; i < D; i++)
    if (!(0.0-epsilon <= xi[i] && xi[i] <= 1.0+epsilon))
      pointIsInElement = false;
    
#if 0  
  tEnd = std::chrono::steady_clock::now();
  durationNewton += std::chrono::duration_cast<std::chrono::duration<double> >(tEnd - tStart).count();
  
  if (nMeasurements%10000 == 0)
  {
    LOG(INFO) << "duration initial approximation: " << durationApproximation / nMeasurements;
    LOG(INFO) << "duration Newton: " << durationNewton / nMeasurements;
    LOG(INFO) << "total: " << (durationNewton+durationApproximation) / nMeasurements;
  }
#endif  
   VLOG(2) << "  " << (pointIsInElement? "inside" : "outside");
  return pointIsInElement;
};

// regular fixed 1D
template<typename BasisFunctionType>
bool FunctionSpaceXi<Mesh::StructuredRegularFixedOfDimension<1>, BasisFunctionType>::
pointIsInElement(Vec3 point, element_no_t elementNo, std::array<double,1> &xi)
{
  const int nDofsPerElement = FunctionSpaceBaseDim<1,BasisFunctionType>::nDofsPerElement();  //=2
  std::array<Vec3, nDofsPerElement> geometryValues;
  
  this->getElementGeometry(elementNo, geometryValues);
  const double elementLength = this->meshWidth() * FunctionSpaceBaseDim<1,BasisFunctionType>::averageNNodesPerElement();
  
  xi[0] = MathUtility::norm<3>(point - geometryValues[0]) / elementLength;
  
  return 0.0 <= xi[0] && xi[0] <= 1.0;
};

// regular fixed 2D
template<typename BasisFunctionType>
bool FunctionSpaceXi<Mesh::StructuredRegularFixedOfDimension<2>, BasisFunctionType>::
pointIsInElement(Vec3 point, element_no_t elementNo, std::array<double,2> &xi)
{
  const int nDofsPerElement = FunctionSpaceBaseDim<2,BasisFunctionType>::nDofsPerElement();
  std::array<Vec3, nDofsPerElement> geometryValues;
  
  this->getElementGeometry(elementNo, geometryValues);
  const double elementLength = this->meshWidth() * FunctionSpaceBaseDim<1,BasisFunctionType>::averageNNodesPerElement();
  
  const double xi1 = (point[0] - geometryValues[0][0]) / elementLength;
  const double xi2 = (point[1] - geometryValues[0][1]) / elementLength;
  
  return (0.0 <= xi1 && xi1 <= 1.0) && (0.0 <= xi2 && xi2 <= 1.0);
};

// regular fixed 3D
template<typename BasisFunctionType>
bool FunctionSpaceXi<Mesh::StructuredRegularFixedOfDimension<3>, BasisFunctionType>::
pointIsInElement(Vec3 point, element_no_t elementNo, std::array<double,3> &xi)
{
  const int nDofsPerElement = FunctionSpaceBaseDim<3,BasisFunctionType>::nDofsPerElement();
  std::array<Vec3, nDofsPerElement> geometryValues;
  
  this->getElementGeometry(elementNo, geometryValues);
  const double elementLength = this->meshWidth() * FunctionSpaceBaseDim<1,BasisFunctionType>::averageNNodesPerElement();
  
  const double xi1 = (point[0] - geometryValues[0][0]) / elementLength;
  const double xi2 = (point[1] - geometryValues[0][1]) / elementLength;
  const double xi3 = (point[2] - geometryValues[0][2]) / elementLength;
  
  return (0.0 <= xi1 && xi1 <= 1.0) && (0.0 <= xi2 && xi2 <= 1.0) && (0.0 <= xi3 && xi3 <= 1.0);
};

// 1D deformable meshes and linear shape function
template<typename MeshType>
bool FunctionSpaceXi<MeshType, BasisFunction::LagrangeOfOrder<1>, Mesh::isDeformableWithDim<1,MeshType>>::
pointIsInElement(Vec3 point, element_no_t elementNo, std::array<double,1> &xi)
{
  //const int nDofsPerElement = FunctionSpaceBaseDim<1,BasisFunction::LagrangeOfOrder<1>>::nDofsPerElement();  //=2
  const int nDofsPerElement = 2;
  std::array<Vec3, nDofsPerElement> geometryValues;
  
  // get needed variables
  this->getElementGeometry(elementNo, geometryValues);
  
  const double xp1 = point[0];
  
  const double x11 = geometryValues[0][0];
  const double x21 = geometryValues[1][0];
  
  // compute analytic solution for xi
  const double xi1 = (x11 - xp1)/(x11 - x21);
  xi[0] = xi1;
  
  return 0.0 <= xi1 && xi1 <= 1.0;
};

// 2D deformable meshes and linear shape function
template<typename MeshType>
bool FunctionSpaceXi<MeshType, BasisFunction::LagrangeOfOrder<1>, Mesh::isDeformableWithDim<2,MeshType>>::
pointIsInElement(Vec3 point, element_no_t elementNo, std::array<double,2> &xi)
{
  //const int nDofsPerElement = FunctionSpaceBaseDim<2,BasisFunction::LagrangeOfOrder<1>>::nDofsPerElement();  //=4
  const int nDofsPerElement = 4;
  std::array<Vec3, nDofsPerElement> geometryValues;
  
  // get needed variables
  this->getElementGeometry(elementNo, geometryValues);
  
  MathUtility::quadrilateralGetPointCoordinates(geometryValues, point, xi);
 
  const double xi1 = xi[0];
  const double xi2 = xi[1];

  return (0.0 <= xi1 && xi1 <= 1.0) && (0.0 <= xi2 && xi2 <= 1.0);
};

// 3D deformable meshes and linear shape function

//! check if the point lies inside the 4-point tetrahedron, if yes, return true and set xi to the value of the point
template<typename MeshType>
bool ComputeXiApproximation<MeshType, BasisFunction::LagrangeOfOrder<1>, Mesh::isDeformableWithDim<3,MeshType>>::
pointIsInTetrahedron(Vec3 point, std::array<Vec3,4> tetrahedron, std::array<bool,3> correctOrientation, std::array<double,3> &xi)
{
  const double xp1 = point[0];
  const double xp2 = point[1];
  const double xp3 = point[2];
  
  const double x11 = tetrahedron[1][0];
  const double x12 = tetrahedron[1][1];
  const double x13 = tetrahedron[1][2];
  
  const double x21 = tetrahedron[2][0];
  const double x22 = tetrahedron[2][1];
  const double x23 = tetrahedron[2][2];
  
  const double x31 = tetrahedron[3][0];
  const double x32 = tetrahedron[3][1];
  const double x33 = tetrahedron[3][2];
  
  const double x41 = tetrahedron[0][0];
  const double x42 = tetrahedron[0][1];
  const double x43 = tetrahedron[0][2];
  
  const double det = (x11 - x41)*(x22 - x42)*(x33 - x43) - (x11 - x41)*(x23 - x43)*(x32 - x42) - (x12 - x42)*(x21 - x41)*(x33 - x43) + (x12 - x42)*(x23 - x43)*(x31 - x41) + (x13 - x43)*(x21 - x41)*(x32 - x42) - (x13 - x43)*(x22 - x42)*(x31 - x41);
  xi[0] = 1./det * ((-x41 + xp1)*((x22 - x42)*(x33 - x43) - (x23 - x43)*(x32 - x42)) + (-x42 + xp2)*(-(x21 - x41)*(x33 - x43) + (x23 - x43)*(x31 - x41)) + (-x43 + xp3)*((x21 - x41)*(x32 - x42) - (x22 - x42)*(x31 - x41)));
  xi[1] = 1./det * ((-x41 + xp1)*(-(x12 - x42)*(x33 - x43) + (x13 - x43)*(x32 - x42)) + (-x42 + xp2)*((x11 - x41)*(x33 - x43) - (x13 - x43)*(x31 - x41)) + (-x43 + xp3)*(-(x11 - x41)*(x32 - x42) + (x12 - x42)*(x31 - x41)));
  xi[2] = 1./det * ((-x41 + xp1)*((x12 - x42)*(x23 - x43) - (x13 - x43)*(x22 - x42)) + (-x42 + xp2)*(-(x11 - x41)*(x23 - x43) + (x13 - x43)*(x21 - x41)) + (-x43 + xp3)*((x11 - x41)*(x22 - x42) - (x12 - x42)*(x21 - x41)));

  const double epsilon = 1e-12;
  bool pointIsInQuarterSpace = (xi[0] >= -epsilon && xi[1] >= -epsilon && xi[2] >= -epsilon);
  
  if (!correctOrientation[0])
    xi[0] = 1. - xi[0];
  if (!correctOrientation[1])
    xi[1] = 1. - xi[1];
  if (!correctOrientation[2])
    xi[2] = 1. - xi[2];
   
  pointIsInQuarterSpace = (xi[0] >= -epsilon && xi[1] >= -epsilon && xi[2] >= -epsilon);
  
  VLOG(3) << "   xi: " << xi << "(in:" << pointIsInQuarterSpace << ")";
  return pointIsInQuarterSpace;
}

// 3D deformable meshes and linear shape function
template<typename MeshType>
bool ComputeXiApproximation<MeshType, BasisFunction::LagrangeOfOrder<1>, Mesh::isDeformableWithDim<3,MeshType>>::
pointIsInElementQuick(Vec3 point, std::array<Vec3, FunctionSpaceBaseDim<3,BasisFunction::LagrangeOfOrder<1>>::nDofsPerElement()> &geometryValues)
{
  VLOG(3) << "pointIsInElementQuick, point " << point << ", element " << geometryValues;

  // this is a heuristic to check whether point is inside the element
  // for every quadrilateral it is checked if the pyramid with the quadrilateral as base and the point as top is oriented such that the point is on the correct side of the element
  // Note that this heuristic is wrong sometimes! Therefore it does not get used anywhere.
  
  const double epsilon = 1e-12;
  
  const Vec3 p0 = (-geometryValues[0]+point);
  const Vec3 p1 = (-geometryValues[1]+point);
  const Vec3 p2 = (-geometryValues[2]+point);
  const Vec3 p3 = (-geometryValues[3]+point);
  const Vec3 p4 = (-geometryValues[4]+point);
  const Vec3 p5 = (-geometryValues[5]+point);
  const Vec3 p6 = (-geometryValues[6]+point);
  const Vec3 p7 = (-geometryValues[7]+point);
  
  if (MathUtility::normSquared<3>(p0) < epsilon || MathUtility::normSquared<3>(p1) < epsilon || MathUtility::normSquared<3>(p2) < epsilon || MathUtility::normSquared<3>(p3) < epsilon
    || MathUtility::normSquared<3>(p4) < epsilon || MathUtility::normSquared<3>(p5) < epsilon || MathUtility::normSquared<3>(p6) < epsilon || MathUtility::normSquared<3>(p7) < epsilon)
  {
    return true;
  }
 
  // bottom [geometryValues[0],geometryValues[1],geometryValues[3],geometryValues[2]]
  const Vec3 temp30 = (-geometryValues[3]+geometryValues[0]);
  const Vec3 temp01 = (-geometryValues[0]+geometryValues[1]);
  const Vec3 temp12 = (-geometryValues[1]+geometryValues[2]);
  const Vec3 temp32 = (-geometryValues[3]+geometryValues[2]);
  const Vec3 temp20 = (-geometryValues[2]+geometryValues[0]);
  
  const bool v0 = MathUtility::dot(MathUtility::cross(temp30, temp01), p0) >= -epsilon;
  const bool v1 = MathUtility::dot(MathUtility::cross(temp01, temp12), p1) >= -epsilon;
  const bool v2 = MathUtility::dot(MathUtility::cross(temp12, temp32), p3) >= -epsilon;
  const bool v3 = MathUtility::dot(MathUtility::cross(temp32, temp20), p2) >= -epsilon;
  
  // top [geometryValues[4],geometryValues[6],geometryValues[7],geometryValues[5]]
  const Vec3 temp74 = (-geometryValues[7]+geometryValues[4]);
  const Vec3 temp46 = (-geometryValues[4]+geometryValues[6]);
  const Vec3 temp65 = (-geometryValues[6]+geometryValues[5]);
  const Vec3 temp75 = (-geometryValues[7]+geometryValues[5]);
  const Vec3 temp54 = (-geometryValues[5]+geometryValues[4]);
  
  const bool v4 = MathUtility::dot(MathUtility::cross(temp74, temp46), p4) >= -epsilon;
  const bool v5 = MathUtility::dot(MathUtility::cross(temp46, temp65), p6) >= -epsilon;
  const bool v6 = MathUtility::dot(MathUtility::cross(temp65, temp75), p7) >= -epsilon;
  const bool v7 = MathUtility::dot(MathUtility::cross(temp75, temp54), p5) >= -epsilon;
  
  // right [geometryValues[1],geometryValues[5],geometryValues[7],geometryValues[3]]
  const Vec3 temp71 = (-geometryValues[7]+geometryValues[1]);
  const Vec3 temp15 = (-geometryValues[1]+geometryValues[5]);
  const Vec3 temp53 = (-geometryValues[5]+geometryValues[3]);
  const Vec3 temp73 = (-geometryValues[7]+geometryValues[3]);
  const Vec3 temp31 = (-geometryValues[3]+geometryValues[1]);
  
  const bool v8 = MathUtility::dot(MathUtility::cross(temp71, temp15), p1) >= -epsilon;
  const bool v9 = MathUtility::dot(MathUtility::cross(temp15, temp53), p5) >= -epsilon;
  const bool v10 = MathUtility::dot(MathUtility::cross(temp53, temp73), p7) >= -epsilon;
  const bool v11 = MathUtility::dot(MathUtility::cross(temp73, temp31), p3) >= -epsilon;
  
  // left [geometryValues[0],geometryValues[2],geometryValues[6],geometryValues[4]]
  const Vec3 temp60 = (-geometryValues[6]+geometryValues[0]);
  const Vec3 temp02 = (-geometryValues[0]+geometryValues[2]);
  const Vec3 temp24 = (-geometryValues[2]+geometryValues[4]);
  const Vec3 temp64 = (-geometryValues[6]+geometryValues[4]);
  const Vec3 temp40 = (-geometryValues[4]+geometryValues[0]);
  
  const bool v12 = MathUtility::dot(MathUtility::cross(temp60, temp02), p0) >= -epsilon;
  const bool v13 = MathUtility::dot(MathUtility::cross(temp02, temp24), p2) >= -epsilon;
  const bool v14 = MathUtility::dot(MathUtility::cross(temp24, temp64), p6) >= -epsilon;
  const bool v15 = MathUtility::dot(MathUtility::cross(temp64, temp40), p4) >= -epsilon;
  
  // front [geometryValues[0],geometryValues[4],geometryValues[5],geometryValues[1]]
  const Vec3 temp50 = (-geometryValues[5]+geometryValues[0]);
  const Vec3 temp04 = (-geometryValues[0]+geometryValues[4]);
  const Vec3 temp41 = (-geometryValues[4]+geometryValues[1]);
  const Vec3 temp51 = (-geometryValues[5]+geometryValues[1]);
  const Vec3 temp10 = (-geometryValues[1]+geometryValues[0]);
  
  const bool v16 = MathUtility::dot(MathUtility::cross(temp50, temp04), p0) >= -epsilon;
  const bool v17 = MathUtility::dot(MathUtility::cross(temp04, temp41), p4) >= -epsilon;
  const bool v18 = MathUtility::dot(MathUtility::cross(temp41, temp51), p5) >= -epsilon;
  const bool v19 = MathUtility::dot(MathUtility::cross(temp51, temp10), p1) >= -epsilon;
  
  // bconst Vec3 tempck [geometryValues[2],geometryValues[3],geometryValues[7],geometryValues[6]]
  const Vec3 temp72 = (-geometryValues[7]+geometryValues[2]);
  const Vec3 temp23 = (-geometryValues[2]+geometryValues[3]);
  const Vec3 temp36 = (-geometryValues[3]+geometryValues[6]);
  const Vec3 temp76 = (-geometryValues[7]+geometryValues[6]);
  const Vec3 temp62 = (-geometryValues[6]+geometryValues[2]);
  
  const bool v20 = MathUtility::dot(MathUtility::cross(temp72, temp23), p2) >= -epsilon;
  const bool v21 = MathUtility::dot(MathUtility::cross(temp23, temp36), p3) >= -epsilon;
  const bool v22 = MathUtility::dot(MathUtility::cross(temp36, temp76), p7) >= -epsilon;
  const bool v23 = MathUtility::dot(MathUtility::cross(temp76, temp62), p6) >= -epsilon;
  
  
  bool isInside = v0 && v1 && v2 && v3 && v4 && v5 && v6 && v7 && v8 && v9 && v10 && v11 && v12 && v13 && v14 
    && v15 && v16 && v17 && v18 && v19 && v20 && v21 && v22 && v23;
    
  VLOG(3) << (isInside? "inside" : "outside");

  if (!isInside)
  {
    VLOG(1) << " isInside: " << v0 << "," <<  v1 << "," <<  v2 << "," <<  v3 << ", " <<  v4 << "," <<  v5 << "," <<  v6 << "," <<  v7 << ", " <<  v8 << "," <<  v9 << "," <<  v10 << "," <<  v11 << ", " <<  v12 << "," <<  v13 << "," <<  v14 
      << "," <<  v15 << ", " <<  v16 << "," <<  v17 << "," <<  v18 << "," <<  v19 << ", " <<  v20 << "," <<  v21 << "," <<  v22 << "," <<  v23;
  }
  else
  {
    VLOG(3) << " isInside: " << v0 << "," <<  v1 << "," <<  v2 << "," <<  v3 << "," <<  v4 << "," <<  v5 << "," <<  v6 << "," <<  v7 << "," <<  v8 << "," <<  v9 << "," <<  v10 << "," <<  v11 << "," <<  v12 << "," <<  v13 << "," <<  v14 
      << "," <<  v15 << "," <<  v16 << "," <<  v17 << "," <<  v18 << "," <<  v19 << "," <<  v20 << "," <<  v21 << "," <<  v22 << "," <<  v23;
  }

  return isInside;
}

template<typename MeshType,typename BasisFunctionType,typename DummyForTraits>
void ComputeXiApproximation<MeshType,BasisFunctionType,DummyForTraits>::
computeApproximateXiForPoint(Vec3 point, element_no_t elementNo, std::array<double,MeshType::dim()> &xi)
{
  xi.fill(0.0);
}

template<typename MeshType>
void ComputeXiApproximation<MeshType, BasisFunction::LagrangeOfOrder<1>, Mesh::isDeformableWithDim<3,MeshType>>::
computeApproximateXiForPoint(Vec3 point, element_no_t elementNo, std::array<double,3> &xi)
{
  // This computes a fast approximation to xi, which can then be refined by a newton scheme
 
  const int nDofsPerElement = FunctionSpaceBaseDim<3,BasisFunction::LagrangeOfOrder<1>>::nDofsPerElement();  //=8
  std::array<Vec3, nDofsPerElement> geometryValues;
  
  this->getElementGeometry(elementNo, geometryValues);
  
  VLOG(3) << "computeApproximateXiForPoint, point " << point << ", element " << elementNo << geometryValues;
  
  xi.fill(0.0);
  std::array<double,3> xiSum{};
  int nSummands = 0;
  
  // p0
  if (pointIsInTetrahedron(point, std::array<Vec3,4>{geometryValues[0], geometryValues[1], geometryValues[2], geometryValues[4]}, std::array<bool,3>{true,true,true}, xiSum))
  {
    xi += xiSum;
    nSummands++;
  }
  else 
  {
    VLOG(3) << "p0 out";
  }
#if 0
  // p1
  if (pointIsInTetrahedron(point, std::array<Vec3,4>{geometryValues[1], geometryValues[0], geometryValues[5], geometryValues[3]}, std::array<bool,3>{false,true,true}, xiSum))
  {
    std::swap(xi[1],xi[2]);
    xi += xiSum;
    nSummands++;
  }
  else 
  {
    VLOG(3) << "p1 out";
  }
  
  //p2
  if (pointIsInTetrahedron(point, std::array<Vec3,4>{geometryValues[2], geometryValues[3], geometryValues[6], geometryValues[0]}, std::array<bool,3>{true,true,false}, xiSum))
  {
    std::swap(xi[1],xi[2]);
    xi += xiSum;
    nSummands++;
  }
  else 
  {
    VLOG(3) << "p2 out";
  }
#endif    
  // p3
  if (pointIsInTetrahedron(point, std::array<Vec3,4>{geometryValues[3], geometryValues[2], geometryValues[1], geometryValues[7]}, std::array<bool,3>{false,false,true}, xiSum))
  {
    xi += xiSum;
    nSummands++;
  }
  else 
  {
    VLOG(3) << "p3 out";
  }
#if   0
  // p4
  if (pointIsInTetrahedron(point, std::array<Vec3,4>{geometryValues[4], geometryValues[5], geometryValues[0], geometryValues[6]}, std::array<bool,3>{true,false,true}, xiSum))
  {
    std::swap(xi[1],xi[2]);
    xi += xiSum;
    nSummands++;
  }
  else 
  {
    VLOG(3) << "p4 out";
  }
#endif  
  // p5
  if (pointIsInTetrahedron(point, std::array<Vec3,4>{geometryValues[5], geometryValues[4], geometryValues[7], geometryValues[1]}, std::array<bool,3>{false,true,false}, xiSum))
  {
    xi += xiSum;
    nSummands++;
  }
  else 
  {
    VLOG(3) << "p5 out";
  }
  
  // p6
  if (pointIsInTetrahedron(point, std::array<Vec3,4>{geometryValues[6], geometryValues[7], geometryValues[4], geometryValues[2]}, std::array<bool,3>{true,false,false}, xiSum))
  {
    xi += xiSum;
    nSummands++;
  }
  else 
  {
    VLOG(3) << "p6 out";
  }
#if 0  
  // p7
  if (pointIsInTetrahedron(point, std::array<Vec3,4>{geometryValues[7], geometryValues[6], geometryValues[3], geometryValues[5]}, std::array<bool,3>{false,false,false}, xiSum))
  {
    std::swap(xi[1],xi[2]);
    xi += xiSum;
    nSummands++;
  }
  else 
  {
    VLOG(3) << "p7 out";
  }
#endif    

  if (nSummands != 0)
    xi /= nSummands;
}

};  // namespace
