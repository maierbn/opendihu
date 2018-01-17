
namespace Mesh
{

template<int D>
MeshD<D>::MeshD(PyObject *specificSettings) : Mesh(specificSettings)
{
}

template<int D>
constexpr int MeshD<D>::dim()
{
  return D;
}

template<int D>
int MeshD<D>::dimension() const
{
  return D;
}

}  // namespace