
namespace Mesh
{

template<int D>
MeshOfDimension<D>::MeshOfDimension(PythonConfig specificSettings) : Mesh(specificSettings)
{
}

template<int D>
constexpr int MeshOfDimension<D>::dim()
{
  return D;
}

template<int D>
int MeshOfDimension<D>::dimension() const
{
  return D;
}

}  // namespace
