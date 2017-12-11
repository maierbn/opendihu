
namespace Mesh
{

template<unsigned long D>
MeshD<D>::MeshD(PyObject *specificSettings) : Mesh(specificSettings)
{
}

template<unsigned long D>
int MeshD<D>::dimension()
{
  return D;
}

}  // namespace