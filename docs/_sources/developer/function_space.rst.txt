FunctionSpace
=============

A function space is a mesh combined with a basis function. This together defines the shape of elements and their nodes. A function space has two template parameters, the mesh and the basis function. 
The following is an example of a 3D quadratic function space.

.. code-block:: c++
  
  FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>, BasisFunction::LagrangeOfOrder<2>>

The function space also knows the location of the nodes in the physical space. Depending on the Mesh class, it stores all node positions (for ``Mesh::StructuredDeformableOfDimension<D>`` and ``Mesh::UnstructuredDeformableOfDimension<D>``) or only the mesh width (for ``Mesh::RegularFixedOfDimension<D>``).

A function space is partioned over multiple processes, on the local process, only the portion of the local subdomain is stored. The case where there is only one process is just a special case of this general concept. Information about the local portion is given by the meshPartition object.

A function space object appears usually as a shared pointer, since multiple object have to know the function space and with a shared pointer it can be provided to every class without copy.

.. code-block:: c++
  
  // define the function space class
  typedef FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<2>, BasisFunction::LagrangeOfOrder<1>> FunctionSpaceType;  // this is an example
  
  // this is how a function space is usually stored, as shared_ptr
  std::shared_ptr<FunctionSpaceType> functionSpace;
  
  // retrieve the mesh partition via
  // functionSpace->meshPartition()

To create a new function space, use the `meshManager` object. This is a singleton class, the mesh manager is available from the `DihuContext` context object. It stores all existent function spaces and takes care that the same function space won't be created twice.

Create a function space from settings (see :doc:`/settings/mesh` for details). If this function space was already created, it returns a pointer to the already created function space. If it was not yet created, it constructs a new object.
This works for all mesh types.

.. code-block:: c++

  // assuming `context_` is the object of type DihuContext and `specificSettings_` is the python settings object of type `PythonConfig`
  // (can be retrived by `specificSettings_ = context_.getPythonConfig();`)
  std::shared_ptr<FunctionSpaceType> functionSpace 
    = context_.meshManager()->functionSpace<FunctionSpaceType>(specificSettings_);
  
Create a function space with a `Mesh::RegularFixedOfDimension<D>` mesh directly, `nElements` is the local number of elements, in the subdomain,
`physicalExtent` is the physical "size" of the subdomain. You can give the function space a name in the first argument, instead of "name". Later the function space can be retrieved from the mesh manager with the `functionSpace<>` method, as above.

.. code-block:: c++

  std::array<element_no_t, D> nElements;
  std::array<double, D> physicalExtent;
  
  std::shared_ptr<FunctionSpaceType> functionSpace 
    = context_.meshManager()->createFunctionSpace<FunctionSpaceType>("name", nElements, physicalExtent);

Create a function space with a `Mesh::StructuredDeformableOfDimension<D>` mesh directly. `nodePositions` is a vector of all local node positions without ghosts, in the local ordering. `nElementsPerCoordinateDirection` is the local number of elements in the subdomain and `nRanksPerCoordinateDirection` is the global number of ranks in each coordinate direction.

.. code-block:: c++

  const int D = 3;   // or 2 or 1
  std::vector<Vec3> nodePositions;
  const std::array<element_no_t,D> nElementsPerCoordinateDirection;
  const std::array<int,D> nRanksPerCoordinateDirection;

  std::shared_ptr<FunctionSpaceType> functionSpace 
    = context_.meshManager()->createFunctionSpace<FunctionSpaceType>("name", nodePositions, nElementsPerCoordinateDirection, nRanksPerCoordinateDirection);
                         

MeshPartition
-------------

A meshPartition object has the following methods:

.. cpp:function:: int nRanks(int coordinateDirection) const
  
  Number of ranks in a coordinate direction.
  
  :param int coordinateDirection: 
  
.. cpp:function:: element_no_t nElementsLocal() const
  
  Number of elements in the current partition.
  
  
.. cpp:function:: global_no_t nElementsGlobal() const
  
  Number of elements in total.
  
  
.. cpp:function:: dof_no_t nDofsLocalWithGhosts() const
  
  Number of dofs in the local partition.
  
  
.. cpp:function:: dof_no_t nDofsLocalWithoutGhosts() const
  
  Number of dofs in the local partition, without ghosts.
  
  
.. cpp:function:: global_no_t nDofsGlobal() const
  
  Number of dofs in total.
  
  
.. cpp:function:: node_no_t nNodesLocalWithGhosts() const
  
  Number of nodes in the local partition.
  
  
.. cpp:function:: node_no_t nNodesLocalWithoutGhosts() const
  
  Number of nodes in the local partition.
  
  
.. cpp:function:: node_no_t nNodesLocalWithGhosts(int coordinateDirection, int partitionIndex = -1) const
  
  Number of nodes in the local partition specified by partitionIndex or the current partition if partitionIndex == -1.
  
  :param int coordinateDirection: 
  :param int partitionIndex = -1: 
  
.. cpp:function:: node_no_t nNodesLocalWithoutGhosts(int coordinateDirection, int partitionIndex = -1) const
  
  Number of nodes in the partition specified by partitionIndex or the current partition if partitionIndex == -1.
  
  :param int coordinateDirection: 
  :param int partitionIndex = -1: 
  
.. cpp:function:: element_no_t nElementsLocal(int coordinateDirection) const
  
  Number of elments in the local partition.
  
  :param int coordinateDirection: 
  
.. cpp:function:: element_no_t nElementsGlobal(int coordinateDirection) const
  
  Number of elments in total.
  
  :param int coordinateDirection: 
  
.. cpp:function:: int beginElementGlobal(int coordinateDirection) const
  
  Global no of first local element.
  
  :param int coordinateDirection: 
  
.. cpp:function:: global_no_t nNodesGlobal() const
  
  Number of nodes in total.
  
  
.. cpp:function:: global_no_t beginNodeGlobalNatural(int coordinateDirection, int partitionIndex = -1) const
  
  Global no of first local node in the partition specified by partitionIndex or the current partition if partitionIndex == -1.
  
  :param int coordinateDirection: 
  :param int partitionIndex = -1: 
  
.. cpp:function:: global_no_t nNodesGlobal(int coordinateDirection) const
  
  Number of nodes in total.
  
  :param int coordinateDirection: 
  
.. cpp:function:: bool hasFullNumberOfNodes(int coordinateDirection, int partitionIndex = -1) const
  
  Get if there are nodes on both boundaries in the given coordinate direction this is the case if the partition touches the right/top/back boundary. Consider the partition specified by partitionIndex or the current partition if partitionIndex == -1.
  
  :param int coordinateDirection: 
  :param int partitionIndex = -1: 
  
.. cpp:function:: const std::vector<element_no_t> &localSizesOnRanks(int coordinateDirection) const
  
  Get a vector with the local sizes on every rank.
  
  :param int coordinateDirection: 
  
.. cpp:function:: ISLocalToGlobalMapping localToGlobalMappingDofs()
  
  Get the local to global mapping for the current partition, for the dof numbering.
  
  
.. cpp:function:: global_no_t getElementNoGlobalNatural(element_no_t elementNoLocal) const
  
  Get the global natural element no for a local element no.
  
  :param element_no_t elementNoLocal: 
  
.. cpp:function:: global_no_t getNodeNoGlobalNatural(std::array<global_no_t,MeshType::dim()> coordinatesGlobal) const
  
  Get the global natural node no for the global coordinates of this node, this can be combined with getCoordinatesGlobal.
  
  :param std\:\:array<global_no_t,MeshType\:\:dim()>: 
  
.. cpp:function:: global_no_t getNodeNoGlobalPetsc(node_no_t nodeNoLocal) const
  
  Get the node no in global petsc ordering from a local node no.
  
  :param node_no_t nodeNoLocal: 
  
.. cpp:function:: void getDofNoGlobalPetsc(const std::vector<dof_no_t> &dofNosLocal, std::vector<PetscInt> &dofNosGlobalPetsc) const
  
  Transfer the local nos in global dof nos, using the PETSc localToGlobal mapping for the dofs.
  
  :param const std\:\:vector<dof_no_t> &dofNosLocal: 
  :param std\:\:vector<PetscInt> &dofNosGlobalPetsc: 
  
.. cpp:function:: global_no_t getDofNoGlobalPetsc(dof_no_t dofNoLocal) const
  
  Get the global petsc dof no for the local no, using the PETSc localToGlobal mapping for the dofs.
  
  :param dof_no_t dofNoLocal: 
  
.. cpp:function:: std::array<global_no_t,MeshType::dim()> getCoordinatesGlobal(node_no_t nodeNoLocal) const
  
  Get the global node coordinates (x,y,z) of the node given by its local node no. This also works for ghost nodes.
  
  
.. cpp:function:: std::array<int,MeshType::dim()> getCoordinatesLocal(node_no_t nodeNoLocal) const
  
  Get the local coordinates for a local node no, also for ghost nodes. With this method and functionSpace->getNodeNo(coordinatesLocal) it is possible to implement a global-to-local mapping.
  
  
.. cpp:function:: std::array<int,MeshType::dim()> getCoordinatesLocal(std::array<global_no_t,MeshType::dim()> coordinatesGlobal, bool &isOnLocalDomain) const
  
  From global natural coordinates compute the local coordinates, set isOnLocalDomain to true if the node with global coordinates is in the local domain.
  
  :param )> getCoordinatesLocal(std\:\:array<global_no_t: 
  :param MeshType\:\:dim()> coordinatesGlobal: 
  :param bool &isOnLocalDomain: 
  
.. cpp:function:: std::array<int,MeshType::dim()> getElementCoordinatesLocal(element_no_t elementNoLocal) const
  
  Get the local coordinates for a local element no.
  
  
.. cpp:function:: element_no_t getElementNoLocal(std::array<int,MeshType::dim()> elementCoordinates) const
  
  Get the local element no. from coordinates.
  
  :param std\:\:array<int,MeshType\:\:dim()>: 
  
.. cpp:function:: node_no_t getNodeNoLocal(global_no_t nodeNoGlobalPetsc) const
  
  Get the local node no for a global petsc node no, does not work for ghost nodes.
  
  :param global_no_t nodeNoGlobalPetsc: 
  
.. cpp:function:: dof_no_t getDofNoLocal(global_no_t dofNoGlobalPetsc) const
  
  Get the local dof no for a global petsc dof no, does not work for ghost nodes.
  
  :param global_no_t dofNoGlobalPetsc: 
  
.. cpp:function:: template <typename T>void extractLocalNodesWithoutGhosts(std::vector<T> &vector, int nComponents=1) const
  
  From a vector of values of global/natural node numbers remove all that are non-local, nComponents consecutive values for each dof are assumed.
  
  :param std\:\:vector<T> &vector: 
  :param int nComponents=1: 
  
.. cpp:function:: template <typename T>void extractLocalDofsWithoutGhosts(std::vector<T> &values) const
  
  From a vector of values of global/natural dofs remove all that are non-local.
  
  :param std\:\:vector<T> &values: 
  
.. cpp:function:: void extractLocalDofsWithoutGhosts(std::vector<double> &values) const
  
  From a vector of values of global/natural dofs remove all that are non-local.
  
  :param std\:\:vector<double> &values: 
  
.. cpp:function:: int convertRankNoToPartitionIndex(int coordinateDirection, int rankNo)
  
  Get the partition index in a given coordinate direction from the rankNo.
  
  :param int coordinateDirection: 
  :param int rankNo: 
  
.. cpp:function:: void output(std::ostream &stream)
  
  Output to stream for debugging.
  
  :param std\:\:ostream &stream: 
  
.. cpp:function:: const std::vector<PetscInt> &dofNosLocal(bool onlyNodalValues=false) const
  
  Get a vector of local dof nos, range [0,nDofsLocalWithoutGhosts] are the dofs without ghost dofs, the whole vector are the dofs with ghost dofs @param onlyNodalValues: if for Hermite only get every second dof such that derivatives are not returned.
  
  :param bool onlyNodalValues=false: 
  
.. cpp:function:: void getDofNosGlobalNatural(std::vector<global_no_t> &dofNosGlobalNatural) const
  
  Get a vector of global natural dof nos of the locally stored non-ghost dofs, needed for setParameters callback function in cellml adapter.
  
  :param std\:\:vector<global_no_t> &dofNosGlobalNatural: 
  
.. cpp:function:: const std::vector<PetscInt> &ghostDofNosGlobalPetsc() const
  
  Get the global dof nos of the ghost dofs in the local partition.
  
  
.. cpp:function:: void initializeDofNosLocalNaturalOrdering(std::shared_ptr<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>> functionSpace)
  
  Initialize the vector dofNosLocalNaturalOrdering\_, this needs the functionSpace and has to be called before dofNosLocalNaturalOrdering() can be used. If the vector is already initialized by a previous call to this method, it has no effect.
  
  :param std\:\:shared_ptr<FunctionSpace\:\:FunctionSpace<MeshType,BasisFunctionType>> functionSpace: 
  
.. cpp:function:: const std::vector<dof_no_t> &dofNosLocalNaturalOrdering() const
  
  Get a vector of local dof nos in local natural ordering, initializeDofNosLocalNaturalOrdering has to be called beforehand.
  
  
.. cpp:function:: bool isNonGhost(node_no_t nodeNoLocal, int &neighbourRankNo) const
  
  Check if the given dof is owned by the own rank, then return true, if not, neighbourRankNo is set to the rank by which the dof is owned.
  
  :param node_no_t nodeNoLocal: 
  :param int &neighbourRankNo: 
  
.. cpp:function:: void getBoundaryElements(Mesh::face_t face, int &neighbourRankNo, std::array<element_no_t,MeshType::dim()> &nBoundaryElements, std::vector<dof_no_t> &dofNos)
  
  Get information about neighbouring rank and boundary elements for specified face, @param neighbourRankNo: the rank of the neighbouring process that shares the face, @param nElements: Size of one-layer mesh that contains boundary elements that touch the neighbouring process.
  
  :param Mesh\:\:face_t face: 
  :param int &neighbourRankNo: 
  :param std\:\:array<element_no_t,MeshType\:\:dim()> &nBoundaryElements: 
  :param std\:\:vector<dof_no_t> &dofNos: 
  
.. cpp:function:: int neighbourRank(Mesh::face_t face)
  
  Get the rank no of the neighbour in direction face, -1 if there is no such neighbour.
  
  :param Mesh\:\:face_t face: 
  
.. cpp:function:: int ownRankPartitioningIndex(int coordinateDirection)
  
  Get the partitioning index in the coordinate direction, i.e. the no. of this rank in this direction, the total number of ranks in each direction can be retrieved by nRanks.
  
  :param int coordinateDirection: 
  


FunctionSpace
-------------

The function space object has the following methods:


.. cpp:function:: static constexpr int nDofsPerElement()
  
  Number of degrees of freedom of this basis.
  
  
.. cpp:function:: static constexpr int nNodesPerElement()
  
  Number of nodes per element.
  
  
.. cpp:function:: static constexpr int nDofsPerNode()
  
  Number of dofs per node.
  
  
.. cpp:function:: static constexpr int averageNDofsPerElement()
  
  If one assigns every dof to an element it is contained in, the number of degrees of freedom per element (not considering boundary elements).
  
  
.. cpp:function:: static constexpr int averageNNodesPerElement()
  
  If one assigns every node to an element it is contained in, the number of nodes per element (not considering boundary elements).
  
  
.. cpp:function:: static double phi(int dofIndex, std::array<double,MeshType::dim()> xi)
  
  Evaluate the basis function corresponding to element-local dof dofIndex at xi, xi lives in [0,1]^D.
  
  :param int dofIndex: 
  :param std\:\:array<double,MeshType\:\:dim()>: 
  
.. cpp:function:: static double dphi_dxi(int dofIndex, int derivativeIdx, std::array<double,MeshType::dim()> xi)
  
  Evaluate the derivative of Phi(xi) w.r.t xi\_i, where i is given by derivativeIdx, i.e. Phi\_{dofIndex,derivativeIdx}(xi).
  
  :param int dofIndex: 
  :param int derivativeIdx: 
  :param std\:\:array<double,MeshType\:\:dim()>: 
  
.. cpp:function:: static std::array<double,MeshType::dim()> gradPhi(int dofIndex, std::array<double,MeshType::dim()> xi)
  
  Evaluate the first derivative of the basis function corresponding to element-local dof dofIndex at xi, interval for xi is [0,1]^D.
  
  :param int dofIndex:
  :param std\:\:array<double,MeshType\:\:dim()>: 
  
.. cpp:function:: void setMeshPartition(std::shared_ptr<Partition::MeshPartition<FunctionSpace<MeshType,BasisFunctionType>,MeshType>> meshPartition)
  
  Set the partition, call this prior to initialize to not initialize the partition from settings but use the given meshPartition.
  
  :param std\:\:shared_ptr<Partition\:\:MeshPartition<FunctionSpace<MeshType,BasisFunctionType>,MeshType>> meshPartition: 
  
.. cpp:function:: std::shared_ptr<Partition::MeshPartition<FunctionSpace<MeshType,BasisFunctionType>,MeshType>> meshPartition() const
  
  Get the partition.
  
  
.. cpp:function:: std::shared_ptr<Partition::MeshPartitionBase> meshPartitionBase()
  
  Get the partition as pointer of type meshPartitionBase, this is in the itnerface in mesh.
  
  
.. cpp:function:: dof_no_t getDofNo(element_no_t elementNoLocal, int dofIndex) const
  
  Return the local dof number of element-local dof dofIndex of element elementNoLocal.
  
  :param element_no_t elementNoLocal: 
  :param int dofIndex: 
  
.. cpp:function:: node_no_t getNodeNo(element_no_t elementNoLocal, int nodeIndex) const
  
  Return the local node number of element-local node nodeIndex of element with local no elementNoLocal.
  
  :param element_no_t elementNoLocal: 
  :param int nodeIndex: 
  
.. cpp:function:: global_no_t getNodeNoGlobalNatural(global_no_t elementNoLocalGlobalNatural, int nodeIndex) const
  
  Return the global/natural node number of element-local node nodeIndex of element with global no elementNoLocalGlobal.
  
  :param global_no_t elementNoLocalGlobalNatural: 
  :param int nodeIndex: 
  
.. cpp:function:: void getNodeDofs(node_no_t nodeGlobalNo, std::vector<dof_no_t> &dofGlobalNos) const
  
  Get all dofs of a specific node, as vector.
  
  :param node_no_t nodeGlobalNo: 
  :param std\:\:vector<dof_no_t> &dofGlobalNos: 
  
.. cpp:function:: void getNodeDofs(node_no_t nodeGlobalNo, std::array<dof_no_t,FunctionSpaceBaseDim<1,BasisFunctionType>::nDofsPerNode()> &dofGlobalNos) const
  
  Get all dofs of a specific node, as array.
  
  :param node_no_t nodeGlobalNo: 
  :param std\:\:array<dof_no_t,FunctionSpaceBaseDim<1,BasisFunctionType>\:\:nDofsPerNode(: 
  
.. cpp:function:: dof_no_t getNodeDofNo(node_no_t nodeGlobalNo, int dofIndex) const
  
  Get the dof no of the specified dof at the node.
  
  :param node_no_t nodeGlobalNo: 
  :param int dofIndex: 
  
.. cpp:function:: node_no_t getNeighbourNodeNoLocal(node_no_t nodeNoLocal, Mesh::face_t direction) const
  
  Get neighbouring node to nodeNoLocal or -1 if there is no such node, nodeNoLocal has to be a non-ghost local node.
  
  :param node_no_t nodeNoLocal: 
  :param Mesh\:\:face_t direction: 
  
.. cpp:function:: node_no_t getNodeNo(std::array<int,MeshType::dim()> coordinateLocal) const
  
  Get node local no from the local coordinate in natural local numbering.
  
  :param std\:\:array<int,MeshType\:\:dim()>: 
  
.. cpp:function:: std::shared_ptr<FieldVariableBaseFunctionSpaceType> fieldVariable(std::string name)
  
  Return a field variable with given name, this is not implemented for structured meshes since there are no extra stored field variables, only for unstructured meshes is it implemented and then stores field variables that were present in parsed exfiles.
  
  :param std\:\:string name: 
  
.. cpp:function:: dof_no_t getDofNoLocal(std::array<global_no_t,MeshType::dim()> coordinatesGlobal, int nodalDofIndex, bool &isOnLocalDomain)
  
  Get the local dof no. for the global coordinates.
  
  :param std\:\:array<global_no_t,MeshType\:\:dim()> coordinatesGlobal: 
  :param int nodalDofIndex: 
  :param bool &isOnLocalDomain: 
  
.. cpp:function:: double meshWidth() const
  
  Get mesh width (=distance between nodes) of the given coordinate direction.
  
  
.. cpp:function:: node_no_t nNodesLocalWithGhosts() const
  
  Return number of nodes including ghost nodes, i.e. these nodes are known locally but some of them are owned by other ranks.
  
  
.. cpp:function:: node_no_t nNodesLocalWithGhosts(int dimension) const
  
  Return number of nodes in specified coordinate direction.
  
  :param int dimension: 
  
.. cpp:function:: node_no_t nNodesLocalWithoutGhosts() const
  
  Return number of nodes that are owned by this partition.
  
  
.. cpp:function:: node_no_t nNodesLocalWithoutGhosts(int dimension) const
  
  Return number of nodes in specified coordinate direction that are owned by this partition.
  
  :param int dimension: 
  
.. cpp:function:: dof_no_t nDofsLocalWithGhosts() const
  
  Return number of dofs.
  
  
.. cpp:function:: dof_no_t nDofsLocalWithoutGhosts() const
  
  Return number of dofs.
  
  
.. cpp:function:: global_no_t nNodesGlobal(int dimension) const
  
  Return number of nodes in specified coordinate direction for the whole global domain.
  
  :param int dimension: 
  
.. cpp:function:: global_no_t nNodesGlobal() const
  
  Return global number of nodes.
  
  
.. cpp:function:: global_no_t nDofsGlobal() const
  
  Return global number of dofs.
  
  
.. cpp:function:: void getNodePositions(std::vector<double> &nodes) const
  
  Fill a vector with the node position entries, nodes will contain consecutively the (x,y,z) values of just all nodes, i.e. for Hermite not the derivatives.
  
  :param std\:\:vector<double> &nodes: 
  
.. cpp:function:: static void getFaceDofs(Mesh::face_t face, std::array<dof_no_t,FunctionSpaceBaseDim<1,BasisFunctionType>::nDofsPerNode()> &dofIndices)
  
  Get all dof indices of a face, note: dimension in FunctionSpaceBaseDim is current-1 (=0), in this case the dofIndices array has exactly so many entries as there are dofs for a node.
  
  :param Mesh\:\:face_t face: 
  :param std\:\:array<dof_no_t,FunctionSpaceBaseDim<1,BasisFunctionType>\:\:nDofsPerNode(: 
  
.. cpp:function:: static int getNeighbourNodeIndex(int nodeIndex, Mesh::face_t face)
  
  Get the neighbouring elemental node index in given direction inside one element or -1 if there is no such node in the element in that direction.
  
  :param int nodeIndex: 
  :param Mesh\:\:face_t face: 
  
.. cpp:function:: std::array<dof_no_t,FunctionSpaceFunction<MeshType,BasisFunctionType>::nNodesPerElement()>getElementNodeNos(element_no_t elementNo) const
  
  Return an array of all node nos. of the element.
  
  
.. cpp:function:: bool findPosition(Vec3 point, element_no_t &elementNo, int &ghostMeshNo, std::array<double,D> &xi, bool startSearchInCurrentElement, double xiTolerance = 1e-4)
  
  Get the element no and the xi value of the point, return true if the point is inside the mesh or false otherwise. Start search at given elementNo ghostMeshNo: -1 means main mesh, 0-5 means ghost Mesh with respecitve Mesh::face\_t.
  
  :param Vec3 point: 
  :param element_no_t &elementNo: 
  :param int &ghostMeshNo: 
  :param std\:\:array<double,D> &xi: 
  :param bool startSearchInCurrentElement: 
  :param double xiTolerance = 1e-4: 
  
.. cpp:function:: bool pointIsInElement(Vec3 point, element_no_t elementNo, std::array<double,D> &xi, double xiTolerance)
  
  Check if the point lies inside the element, if yes, return true and set xi to the value of the point, defined in 11\_function\_space\_xi.h.
  
  :param Vec3 point: 
  :param element_no_t elementNo: 
  :param std\:\:array<double,D> &xi: 
  :param double xiTolerance: 
  
.. cpp:function:: void setGhostMesh(Mesh::face_t face, const std::shared_ptr<FunctionSpace<MeshType,BasisFunctionType>> ghostMesh)
  
  Store a ghost mesh which is a neighouring mesh with only one layer of elements, this will be used by pointIsInElement and findPosition.
  
  :param Mesh\:\:face_t face: 
  :param const std\:\:shared_ptr<FunctionSpace<MeshType,BasisFunctionType>> ghostMesh: 
  
.. cpp:function:: bool findPosition(Vec3 point, element_no_t &elementNo, int &ghostMeshNo, std::array<double,MeshType::dim()> &xi, bool startSearchInCurrentElement, double xiTolerance = 1e-4)
  
  Get the element no and the xi value of the point, return true if the point is inside the mesh or false otherwise. Start search at given elementNo ghostMeshNo: -1 means main mesh, 0-5 means ghost Mesh with respecitve Mesh::face\_t.
  
  :param Vec3 point: 
  :param element_no_t &elementNo: 
  :param int &ghostMeshNo: 
  :param std\:\:array<double,MeshType\:\:dim()> &xi: 
  :param bool startSearchInCurrentElement: 
  :param double xiTolerance = 1e-4: 
  
.. cpp:function:: bool checkNeighbouringElements(const Vec3 &point, element_no_t &elementNo, int &ghostMeshNo, std::array<double,MeshType::dim()> &xi)
  
  Check if the point is in a neighbouring element to elementNo on ghostMeshNo (-1=main mesh, 0-5=ghost mesh on respective face, 0=face0Minus, 1=face0Plus, etc.), return true if the element was found amoung the neighbours set elementNo, ghostMeshNo and xi appropriately.
  
  :param const Vec3 &point: 
  :param element_no_t &elementNo: 
  :param int &ghostMeshNo: 
  :param std\:\:array<double,MeshType\:\:dim()>: 
  
.. cpp:function:: std::shared_ptr<FieldVariable::FieldVariableBaseFunctionSpace<FunctionSpace<MeshType,BasisFunctionType>>> createFieldVariable(std::string name, std::vector<std::string> componentNames)
  
  Create a non-geometry field field variable with no values being set, with given component names.
  
  :param std\:\:string name: 
  :param std\:\:vector<std\:\:string> componentNames: 
  
.. cpp:function:: std::shared_ptr<FieldVariable::FieldVariableBaseFunctionSpace<FunctionSpace<MeshType,BasisFunctionType>>> createFieldVariable(std::string name, int nComponents=1)
  
  Create a non-geometry field field variable with no values being set, with given number of components, the component names will be the numbers.
  
  :param std\:\:string name: 
  :param int nComponents=1: 
  
.. cpp:function:: template <int nComponents>std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace<MeshType,BasisFunctionType>,nComponents>> createFieldVariable(std::string name)
  
  Create a non-geometry field field variable with no values being set, with given number of components, the component names will be the numbers.
  
  :param std\:\:string name: 
  
.. cpp:function:: template <int nComponents>std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace<MeshType,BasisFunctionType>,nComponents>> createFieldVariable(std::string name, std::vector<std::string> componentNames)
  
  Create a non-geometry field field variable with no values being set, with given number of components and component names.
  
  :param std\:\:string name: 
  :param std\:\:vector<std\:\:string> componentNames: 
  
.. cpp:function:: int getNumberScaleFactors(element_no_t elementGlobalNo)
  
  Get the number of scale factors that are stored for the global element no.
  
  :param element_no_t elementGlobalNo: 
  
.. cpp:function:: std::array<std::array<double,MeshType::dim()>,FunctionSpaceFunction<MeshType,BasisFunctionType>::nDofsPerElement()>getGradPhi(std::array<double,MeshType::dim()> xi) const
  
  Return an array of the gradients of all nodal basis functions, evaluated at xi.
  
  :param )>: 
  :param FunctionSpaceFunction<MeshType,BasisFunctionType>\:\:nDofsPerElement()>getGradPhi(std\:\:array<double: 
  :param MeshType\:\:dim()>: 
  
.. cpp:function:: template <int nComponents>std::array<double,nComponents> interpolateValueInElement(std::array<std::array<double,nComponents>,FunctionSpaceFunction<MeshType,BasisFunctionType>::nDofsPerElement()> &elementalDofValues,std::array<double,MeshType::dim()> xi) const
  
  Interpolate the nComponents values within an element at the given xi position using the basis functions.
  
  :param std\:\:array<std\:\:array<double,nComponents>,FunctionSpaceFunction<MeshType,BasisFunctionType>\:\:nDofsPerElement()> &elementalDofValues: 
  :param std\:\:array<double,MeshType\:\:dim()>: 
  
.. cpp:function:: double interpolateValueInElement(std::array<double,FunctionSpaceFunction<MeshType,BasisFunctionType>::nDofsPerElement()> &elementalDofValues,std::array<double,MeshType::dim()> xi) const
  
  Interpolate the value within an element at the given xi position using the basis functions.
  
  :param std\:\:array<double,FunctionSpaceFunction<MeshType,BasisFunctionType>\:\:nDofsPerElement()> &elementalDofValues: 
  :param std\:\:array<double,MeshType\:\:dim()>: 
  
.. cpp:function:: std::array<double,MeshType::dim()> interpolateGradientInElement(std::array<double,FunctionSpaceFunction<MeshType,BasisFunctionType>::nDofsPerElement()> &elementalDofValues,Tensor2<MeshType::dim()> inverseJacobianParameterSpace, std::array<double,MeshType::dim()> xi) const
  
  Interpolate the gradient of a scalar field within an element at the given xi position using the basis functions the inverseJacobianParameterSpace can be computed by getInverseJacobian.
  
  :param )> interpolateGradientInElement(std\:\:array<double: 
  :param FunctionSpaceFunction<MeshType,BasisFunctionType>\:\:nDofsPerElement()> &elementalDofValues: 
  :param Tensor2<MeshType\:\:dim()> inverseJacobianParameterSpace: 
  :param std\:\:array<double,MeshType\:\:dim()>: 
  
.. cpp:function:: Vec3 getNormal(Mesh::face_t face, std::array<Vec3,FunctionSpaceFunction<MeshType,BasisFunctionType>::nDofsPerElement()> geometryValues, std::array<double,MeshType::dim()> xi)
  
  Compute the normal in world space, normal to face at xi, use the given geometry values, that can by obtained by fieldVariable->getElementValues(elementNo, geometryValues) or mesh->getElementGeometry(elementNo, geometryValues).
  
  :param Mesh\:\:face_t face: 
  :param std\:\:array<Vec3,FunctionSpaceFunction<MeshType,BasisFunctionType>\:\:nDofsPerElement()> geometryValues: 
  :param std\:\:array<double,MeshType\:\:dim()>: 
  
.. cpp:function:: Vec3 getNormal(Mesh::face_t face, element_no_t elementNoLocal, std::array<double,MeshType::dim()> xi)
  
  Compute the normal in world space, normal to face at xi.
  
  :param Mesh\:\:face_t face: 
  :param element_no_t elementNoLocal: 
  :param std\:\:array<double,MeshType\:\:dim()>: 
  
.. cpp:function:: Tensor2<MeshType::dim()> getInverseJacobian(std::array<Vec3,FunctionSpaceFunction<MeshType,BasisFunctionType>::nDofsPerElement()> &geometryValues, element_no_t elementNo, std::array<double,MeshType::dim()> xi)
  
  Compute the inverseJacobian that is needed to transform a gradient vector from parameter space to world space, for an element at a xi position. This version of the method needs the values of the geometry field, if the jacobian is needed at multiple positions in the same element, these values can be retrieved once and used for all computations of the jacobians. There is also the convienience method which does not need the geometryValues but gets them itself. The following properties of the jacobian hold: jacobianParameterSpace[columnIdx][rowIdx] = dX\_rowIdx/dxi\_columnIdx inverseJacobianParameterSpace[columnIdx][rowIdx] = dxi\_rowIdx/dX\_columnIdx because of inverse function theorem.
  
  :param )> getInverseJacobian(std\:\:array<Vec3: 
  :param FunctionSpaceFunction<MeshType,BasisFunctionType>\:\:nDofsPerElement()> &geometryValues: 
  :param element_no_t elementNo: 
  :param std\:\:array<double,MeshType\:\:dim()>: 
  
.. cpp:function:: Tensor2<MeshType::dim()> getInverseJacobian(element_no_t elementNo, std::array<double,MeshType::dim()> xi)
  
  Compute the inverseJacobian that is needed to transform a gradient vector from parameter space to world space, for an element at a xi position. The following properties of the jacobian hold: jacobianParameterSpace[columnIdx][rowIdx] = dX\_rowIdx/dxi\_columnIdx inverseJacobianParameterSpace[columnIdx][rowIdx] = dxi\_rowIdx/dX\_columnIdx because of inverse function theorem.
  
  :param )> getInverseJacobian(element_no_t elementNo: 
  :param std\:\:array<double,MeshType\:\:dim()>: 
  
.. cpp:function:: bool pointIsInElement(Vec3 point, element_no_t elementNo, std::array<double,1> &xi, double xiTolerance = 1e-4)
  
  Check if the point lies inside the element, if yes, return true and set xi to the value of the point.
  
  :param Vec3 point: 
  :param element_no_t elementNo: 
  :param std\:\:array<double,1> &xi: 
  :param double xiTolerance = 1e-4: 
  
.. cpp:function:: std::array<dof_no_t,FunctionSpaceFunction<MeshType,BasisFunctionType>::nDofsPerElement()>getElementDofNosLocal(element_no_t elementNo) const
  
  Return an array of all dof nos. of the element, including ghost dofs (local dof nos).
  
  
.. cpp:function:: void getElementDofNosLocalWithoutGhosts(element_no_t elementNo, std::vector<dof_no_t> &dofNosLocal) const
  
  Fill a vector of all local dof nos. of the element, without ghost dofs.
  
  :param element_no_t elementNo: 
  :param std\:\:vector<dof_no_t> &dofNosLocal: 
  

.. cpp:function:: Vec3 getGeometry(node_no_t dofGlobalNo) const
  
  Return the geometry field entry (node position for Lagrange elements) of a specific dof.
  
  :param node_no_t dofGlobalNo: 
  
.. cpp:function:: void getElementGeometry(element_no_t elementNoLocal, std::array<Vec3, FunctionSpaceBaseDim<MeshType::dim(),BasisFunctionType>::nDofsPerElement()> &values)
  
  Get all geometry entries for an element.
  
  :param element_no_t elementNoLocal: 
  :param std\:\:array<Vec3,FunctionSpaceBaseDim<MeshType\:\:dim(),BasisFunctionType>\:\:nDofsPerElement(: 
  
.. cpp:function:: void extractSurfaceGeometry(const std::array<Vec3, FunctionSpaceBaseDim<MeshType::dim(),BasisFunctionType>::nDofsPerElement()> &geometryVolume, Mesh::face_t face,std::array<Vec3, FunctionSpaceBaseDim<MeshType::dim()-1,BasisFunctionType>::nNodesPerElement()> &geometrySurface)
  
  From the function space geometry, extract geometry data for a surface with has one lower dimensionality, only the nodal dofs are extracted, also for Hermite.
  
  :param const std\:\:array<Vec3,FunctionSpaceBaseDim<MeshType\:\:dim(),BasisFunctionType>\:\:nDofsPerElement()> &geometryVolume: 
  :param Mesh\:\:face_t face: 
  :param std\:\:array<Vec3,FunctionSpaceBaseDim<MeshType\:\:dim()-1,BasisFunctionType>\:\:nNodesPerElement(: 
  
.. cpp:function:: GeometryFieldType &geometryField()
  
  Return the internal geometry field variable.
  
  


