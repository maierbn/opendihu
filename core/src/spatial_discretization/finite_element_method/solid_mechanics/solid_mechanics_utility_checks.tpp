#include "spatial_discretization/finite_element_method/solid_mechanics/solid_mechanics_utility.h"

#include "easylogging++.h"
#include "utility/math_utility.h"
#include "control/types.h"

#include <cmath>
#include <array>

namespace SpatialDiscretization
{

template<int D>
void checkSymmetry(const Tensor2<D> &rightCauchyGreen, std::string name)
{
  if (VLOG_IS_ON(1))
  {
   
    bool isSymmetric = true;
    const double errorTolerance = 1e-14;
    for (int a=0; a<D; a++)
    {
      for (int b=0; b<D; b++)
      {
        if (fabs(rightCauchyGreen[a][b] - rightCauchyGreen[b][a]) > errorTolerance)
        {
          LOG(ERROR) << name << "["<<a<<"]["<<b<<"] != " << name << "["<<b<<"]["<<a<<"] ("<<rightCauchyGreen[a][b]<<" != "<<rightCauchyGreen[b][a]<<") - symmetry violated"<<std::endl;
          isSymmetric = false;
        }
      }
    }
    if (isSymmetric)
      VLOG(2) << name << " is symmetric";
  }
}

template<int D>
void checkInverseIsCorrect(const Tensor2<D> &rightCauchyGreen, Tensor2<D> &inverseRightCauchyGreen, std::string name)
{
 
  if (VLOG_IS_ON(1))
  {
    bool inverseCorrect = true;
    double errorTolerance = 1e-12;
    for (int a=0; a<D; a++)
    {
      for (int b=0; b<D; b++)
      {
        double matrixProduct = 0.0;
        for (int k=0; k<D; k++)
        {
          matrixProduct += rightCauchyGreen[k][a] * inverseRightCauchyGreen[b][k];
        }
        double delta_ab = (a == b? 1.0 : 0.0);
        
        if (fabs(delta_ab - matrixProduct) > errorTolerance)
        {
          LOG(ERROR) << name << " or " << name <<"^{-1} is wrong: " << matrixProduct << " should be " << delta_ab;
          inverseCorrect = false;
        }
      }
    }
    
    if (inverseCorrect)
    {
      VLOG(2) << name << " corresponds to " << name <<"^{-1}";
    }
    else
    {
      LOG(DEBUG) << name << ": " << rightCauchyGreen << std::endl << ", inv: " << inverseRightCauchyGreen << std::endl;
    }
  }
}

};  // namespace