#include "mesh/structured_deformable.h"

namespace Mesh
{

template<int D>
StructuredDeformableOfDimension<D>::
StructuredDeformableOfDimension(PythonConfig specificSettings) :
  Structured<D>(specificSettings), hasTriangleCorners_(false)
{
  if (specificSettings.hasKey("hasTriangleCorners"))
  {
    if (D == 3)
    {
#ifdef USE_VECTORIZED_FE_MATRIX_ASSEMBLY
      LOG(WARNING) << "Option \"hasTriangleCorners\" was specified when USE_VECTORIZED_FE_MATRIX_ASSEMBLY is set. This is not implemented.";
#else
      hasTriangleCorners_ = specificSettings.getOptionBool("hasTriangleCorners", false);
#endif
    }
    else
    {
      LOG(WARNING) << "Option \"hasTriangleCorners\" was specified for a mesh with D=" << D << ". It is only possible for D=3.";
    }
  }
}

template<int D>
bool StructuredDeformableOfDimension<D>::
hasTriangleCorners()
{
  return hasTriangleCorners_;
}

}  // namespace
