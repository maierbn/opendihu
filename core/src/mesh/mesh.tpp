
namespace Mesh
{

template<unsigned long D>
MeshD<D>::MeshD(PyObject *specificSettings) : Mesh(specificSettings)
{
}

template<unsigned long D>
constexpr int MeshD<D>::dim()
{
  return D;
}

template<unsigned long D>
int MeshD<D>::dimension() const
{
  return D;
}

}  // namespace