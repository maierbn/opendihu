FieldVariable
=============

Field variables are the discretized scalar or vector valued functions over the domain. They depend on the function space and have a fixed number of components.

.. code-block:: c++
  
  FieldVariable::FieldVariable<FunctionSpaceType,nComponents>


A field variable is a vector with as many entries as there are unknowns in the function space. Each entry can have multiple components. 
For example, a geometry field variable, that stores the geometric position of each node of the mesh has 3 components for the `x`, `y` and `z` values.

The field variable itself has a name, that is used when it is written to output files, e.g. to identify the variable in Paraview, and for debugging output.
Also the components have names.

Internally, field variables use a PartitionedPetscVec as data container which is a wrapper for a Petsc Vec.

Field variables can be generated from the function space on which they are defined:

.. code-block:: c++
  
  // define a field variable with 1 component (scalar)
  // and name "flowPotential"
  fieldVariable_ = functionSpace_->template createFieldVariable<1>("flowPotential");

  // vector-valued field variable with 3 components x,y,z
  std::vector<std::string> componentNames({"x","y","z"});
  fieldVariable2_ = functionSpace_->template createFieldVariable<3>("Δu", componentNames);
  

General methods
------------------

.. cpp:function:: void outputHeaderExelem(std::ostream &file, element_no_t currentElementGlobalNo, int fieldVariableNo=-1)
  
  Write a exelem file header to a stream, for a particular element, fieldVariableNo is the field index x) in the exelem file header. For parallel program execution this writes headers for the local exelem files on every rank.
  
  :param std\:\:ostream &file: 
  :param element_no_t currentElementGlobalNo: 
  :param int fieldVariableNo=-1: 
  
.. cpp:function:: void outputHeaderExnode(std::ostream &file, node_no_t currentNodeGlobalNo, int &valueIndex, int fieldVariableNo=-1)
  
  Write a exelem file header to a stream, for a particular node, for parallel program execution this writes headers for the local exnodes files on every rank.
  
  :param std\:\:ostream &file: 
  :param node_no_t currentNodeGlobalNo: 
  :param int &valueIndex: 
  :param int fieldVariableNo=-1: 
  
.. cpp:function:: bool haveSameExfileRepresentation(element_no_t element1, element_no_t element2)
  
  Tell if 2 elements have the same exfile representation, i.e. same number of versions.
  
  :param element_no_t element1: 
  :param element_no_t element2: 
  
.. cpp:function:: Vec &valuesLocal(int componentNo = 0)
  
  Get the internal PETSc vector values, the local vector for the specified component.
  
  :param int componentNo = 0: 
  
.. cpp:function:: Vec &valuesGlobal(int componentNo)
  
  Get the internal PETSc vector values, the global vector for the specified component.
  
  :param int componentNo: 
  
.. cpp:function:: Vec &valuesGlobal()
  
  If the vector has multiple components, return a nested Vec of the global vector, else return the global vector.
  
  
.. cpp:function:: Vec &getValuesContiguous()
  
  Fill a contiguous vector with all components after each other, "struct of array"-type data layout. after manipulation of the vector has finished one has to call restoreValuesContiguous.
  
  
.. cpp:function:: void restoreValuesContiguous()
  
  Copy the values back from a contiguous representation where all components are in one vector to the standard internal format of PartitionedPetscVec where there is one local vector with ghosts for each component. this has to be called.
  
  
.. cpp:function:: void output(std::ostream &stream) const
  
  Output string representation to stream for debugging.
  
  :param std\:\:ostream &stream: 
  
.. cpp:function:: std::shared_ptr<FunctionSpaceType> functionSpace()
  
  Return the functionSpace of this field variable.
  
  
.. cpp:function:: std::string name() const
  
  Get the name of the field variable.
  
  
.. cpp:function:: bool isGeometryField() const
  
  If the field has the flag "geometry field", i.e. in the exelem file its type was specified as "coordinate".
  
  
.. cpp:function:: void checkNansInfs(int componentNo = 0) const
  
  Check if there are NaNs or high values in the current variable, if yes output a warning.
  
  :param int componentNo = 0: 
  
.. cpp:function:: dof_no_t nDofsLocalWithoutGhosts() const
  
  Get the number of dofs.
  
  
.. cpp:function:: dof_no_t nDofsGlobal() const
  
  Get the number of global dofs.
  
  
.. cpp:function:: const std::array<std::string,nComponentsValue> &componentNames() const
  
  Get the component names.
  
  
.. cpp:function:: virtual std::shared_ptr<Component<FunctionSpaceType,nComponentsValue>> component(int componentNo) = 0
  
  Return the component by index.
  
  :param int componentNo: 
  
.. cpp:function:: const std::string componentName(int componentNo) const
  
  Get the component Name.
  
  :param int componentNo: 
  
.. cpp:function:: static constexpr int nComponents()
  
  Get the number of components.
  
  
.. cpp:function:: virtual int getNComponents() const
  
  Get the number of components.
  
  
.. cpp:function:: void computeGradientField(std::shared_ptr<FieldVariable<FunctionSpaceType, FunctionSpaceType::dim()>> gradientField,std::shared_ptr<FieldVariable<FunctionSpaceType,1>> jacobianConditionNumber = nullptr)
  
  Fill the gradient field with the gradient values in world coordinates of this field variable. This is only possible for scalar fields.
  
  :param std\:\:shared_ptr<FieldVariable<FunctionSpaceType,FunctionSpaceType\:\:dim()>> gradientField: 
  :param std\:\:shared_ptr<FieldVariable<FunctionSpaceType,1>> jacobianConditionNumber = nullptr: 
  
.. cpp:function:: void startGhostManipulation()
  
  This has to be called before the vector is manipulated (i.e. VecSetValues or vecZeroEntries is called), to ensure that the current state of the vector is fetched from the global vector.
  
  
.. cpp:function:: void zeroGhostBuffer()
  
  Zero all values in the local ghost buffer. Needed if between startGhostManipulation() and finishGhostManipulation() only some ghost will be reassigned. To prevent that the "old" ghost values that were present in the local ghost values buffer get again added to the real values which actually did not change.
  
  
.. cpp:function:: void finishGhostManipulation()
  
  This has to be called after the vector is manipulated (i.e. VecSetValues or vecZeroEntries is called), to ensure that operations on different partitions are merged by Petsc It sums up the values in the ghost buffer and the actual nodal value.
  
  
.. cpp:function:: void setRepresentationGlobal()
  
  Set the internal representation to be global, i.e. using the global vectors, if it was local, ghost buffer entries are discarded (use finishGhostManipulation to consider ghost dofs).
  
  
.. cpp:function:: void setRepresentationLocal()
  
  Set the internal representation to be local, i.e. using the local vectors, ghost buffer is not filled (use startGhostManipulation to consider ghost dofs).
  
  
.. cpp:function:: void setRepresentationContiguous()
  
  Set the internal representation to be contiguous, i.e. using the contiguous vectors for a specific component, get all values.
  
  :param );void getValuesWithGhosts(int componentNo: 
  :param std\:\:vector<double> &values: 
  :param bool onlyNodalValues=false: if this is true, for Hermite only the non-derivative values are retrieved.
  
Getters
------------
  
.. cpp:function:: void getValuesWithoutGhosts(int componentNo, std::vector<double> &values, bool onlyNodalValues=false) const
  
  For a specific component, get all values
  
  :param int componentNo: 
  :param std\:\:vector<double> &values: 
  :param bool onlyNodalValues=false: if this is true, for Hermite only the non-derivative values are retrieved.
  
.. cpp:function:: void getValuesWithGhosts(std::vector<std::array<double,nComponents>> &values, bool onlyNodalValues=false) const
  
  Get all values
  
  :param std\:\:vector<std\:\:array<double,nComponents>> &values: 
  :param bool onlyNodalValues=false: if this is true, for Hermite only the non-derivative values are retrieved.
  
.. cpp:function:: void getValuesWithoutGhosts(std::vector<std::array<double,nComponents>> &values, bool onlyNodalValues=false) const
  
  Get all values
  
  :param std\:\:vector<std\:\:array<double,nComponents>> &values: 
  :param bool onlyNodalValues=false: if this is true, for Hermite only the non-derivative values are retrieved.
  
.. cpp:function:: void getValuesWithoutGhosts(std::array<std::vector<double>,nComponents> &values, bool onlyNodalValues=false) const
  
  Get all values
  
  :param std\:\:array<std\:\:vector<double>,nComponents> &values: 
  :param bool onlyNodalValues=false: if this is true, for Hermite only the non-derivative values are retrieved.
  
.. cpp:function:: template<int N>void getValues(int componentNo, std::array<dof_no_t,N> dofLocalNo, std::array<double,N> &values) const
  
  For a specific component, get values from their local dof no.s, as array, therefore templated by the number of elements, N, to retrieve.
  
  :param int componentNo: 
  :param std\:\:array<dof_no_t,N> dofLocalNo: 
  :param std\:\:array<double,N> &values: 
  
.. cpp:function:: void getValues(int componentNo, const std::vector<dof_no_t> &dofLocalNo, std::vector<double> &values) const
  
  For a specific component, get values from their local dof no.s, as vector.
  
  :param int componentNo: 
  :param const std\:\:vector<dof_no_t> &dofLocalNo: 
  :param std\:\:vector<double> &values: 
  
.. cpp:function:: void getValues(const std::vector<dof_no_t> &dofLocalNo, std::vector<double> &values) const
  
  Get values for all components, from their local dof no.s, as contiguous vector in order [comp0, comp0, comp0, ..., comp1, comp1, ...].
  
  :param const std\:\:vector<dof_no_t> &dofLocalNo: 
  :param std\:\:vector<double> &values: 
  
.. cpp:function:: template<int N>void getValues(std::array<dof_no_t,N> dofLocalNo, std::array<std::array<double,nComponents>,N> &values) const
  
  Get values from their local dof no.s for all components.
  
  :param std\:\:array<dof_no_t,N> dofLocalNo: 
  :param std\:\:array<std\:\:array<double,nComponents>,N> &values: 
  
.. cpp:function:: void getValues(std::vector<dof_no_t> dofLocalNo, std::vector<std::array<double,nComponents>> &values) const
  
  Get values from their local dof no.s for all components.
  
  :param std\:\:vector<dof_no_t> dofLocalNo: 
  :param std\:\:vector<std\:\:array<double,nComponents>> &values: 
  
.. cpp:function:: void getElementValues(int componentNo, element_no_t elementNoLocal, std::array<double,FunctionSpaceType::nDofsPerElement()> &values) const
  
  For a specific component, get the values corresponding to all element-local dofs.
  
  :param int componentNo: 
  :param element_no_t elementNoLocal: 
  :param std\:\:array<double,FunctionSpaceType\:\:nDofsPerElement(: 
  
.. cpp:function:: void getElementValues(element_no_t elementNoLocal, std::array<std::array<double,nComponents>,FunctionSpaceType::nDofsPerElement()> &values) const
  
  Get the values corresponding to all element-local dofs for all components.
  
  :param element_no_t elementNoLocal: 
  :param std\:\:array<std\:\:array<double,nComponents>,FunctionSpaceType\:\:nDofsPerElement(: 
  
.. cpp:function:: double getValue(int componentNo, node_no_t dofLocalNo) const
  
  For a specific component, get a single value from local dof no.
  
  :param int componentNo: 
  :param node_no_t dofLocalNo: 
  
.. cpp:function:: getValues(int componentNo, int nValues, const dof_no_t *dofLocalNo, std::vector<double> &values) const
  
  Get values from their local dof no.s for a specific component.
  
  :param int componentNo:
  :param int nValues: The number of values to get.
  :param const dof_no_t \*dofLocalNo:
  :param std\:\:vector<double> &values: The resulting values will be appended to this vector.
  
.. cpp:function:: std::array<double,nComponents> getValue(node_no_t dofLocalNo) const
  
  Get a single value from local dof no. for all components.
  
  :param node_no_t dofLocalNo: 
  
.. cpp:function:: void extractComponentCopy(int componentNo, std::shared_ptr<FieldVariable<FunctionSpaceType,1>> extractedFieldVariable)
  
  Extract the specified component from the field variable (by copying it) and store it in the given field variable (which already has the data allocated).
  
  :param int componentNo: 
  :param std\:\:shared_ptr<FieldVariable<FunctionSpaceType,1>> extractedFieldVariable: 
  
.. cpp:function:: void extractComponentShared(int componentNo, std::shared_ptr<FieldVariable<FunctionSpaceType,1>> extractedFieldVariable)
  
  Extract the specified component from the field variable by using the raw data array in the given field variable. Afterwards this field variable is invalid and can only be used again after restoreExtractedComponent has been called.
  
  :param int componentNo: 
  :param std\:\:shared_ptr<FieldVariable<FunctionSpaceType,1>> extractedFieldVariable: 
  
.. cpp:function:: template<int nComponents2>void restoreExtractedComponent(std::shared_ptr<PartitionedPetscVec<FunctionSpaceType,nComponents2>> extractedVec)
  
  Restore the extracted raw array to petsc and make the field variable usable again.
  
  :param std\:\:shared_ptr<PartitionedPetscVec<FunctionSpaceType,nComponents2>> extractedVec: 
  
.. cpp:function:: void getElementValues(element_no_t elementNoLocal, std::array<double,FunctionSpaceType::nDofsPerElement()> &values) const
  
  Only for scalar field variables:  get the values corresponding to all element-local dofs for all components.
  
  :param element_no_t elementNoLocal: 
  :param std\:\:array<double,FunctionSpaceType\:\:nDofsPerElement(: 
  
.. cpp:function:: double getValue(node_no_t dofLocalNo) const
  
  Only for scalar field variables:  get a single value from local dof no. for all components.
  
  :param node_no_t dofLocalNo: 
  
.. cpp:function:: void getValuesWithGhosts(std::vector<double> &values, bool onlyNodalValues=false) const
  
  Only for scalar field variables:  get all stored local values.
  
  :param std\:\:vector<double> &values: 
  :param bool onlyNodalValues=false: 
  
.. cpp:function:: void getValuesWithoutGhosts(std::vector<double> &values, bool onlyNodalValues=false) const
  
  Only for scalar field variables:  get all stored local values.
  
  :param std\:\:vector<double> &values: 
  :param bool onlyNodalValues=false: 
  
Setters
------------

.. cpp:function:: void setValues(int componentNo, Vec petscVector)
  
  Set the values for the given component from a petsc Vec.
  
  :param int componentNo: 
  :param Vec petscVector: 
  
.. cpp:function:: void setValues(int componentNo, std::shared_ptr<FieldVariable<FunctionSpaceType,1>> fieldVariable)
  
  Set the values for the given component from the other field variable.
  
  :param int componentNo: 
  :param std\:\:shared_ptr<FieldVariable<FunctionSpaceType,1>> fieldVariable: 
  
.. cpp:function:: void setValues(int componentNo, const std::vector<dof_no_t> &dofNosLocal, const std::vector<double> &values, InsertMode petscInsertMode=INSERT_VALUES)
  
  Set values for a given component for given dofs.
  
  :param int componentNo: 
  :param const std\:\:vector<dof_no_t> &dofNosLocal: 
  :param const std\:\:vector<double> &values: 
  :param InsertMode petscInsertMode=INSERT_VALUES: 
  
.. cpp:function:: template<int N>void setValues(int componentNo, const std::array<dof_no_t,N> &dofNosLocal, const std::array<double,N> &values, InsertMode petscInsertMode=INSERT_VALUES)
  
  Set values for a given component for given dofs.
  
  :param int componentNo: 
  :param const std\:\:array<dof_no_t,N> &dofNosLocal: 
  :param const std\:\:array<double,N> &values: 
  :param InsertMode petscInsertMode=INSERT_VALUES: 
  
.. cpp:function:: void setValues(int componentNo, int nValues, const dof_no_t *dofNosLocal, const double *values, InsertMode petscInsertMode=INSERT_VALUES)
  
  Set values for a given component for given dofs, using raw pointers. This directly mirors the VecSetValues function of Petsc.
  
  :param int componentNo:
  :param int nValues: The number of values to set
  :param const dof_no_t \*dofNosLocal:
  :param const double \*values:
  :param InsertMode petscInsertMode=INSERT_VALUES: 
  
.. cpp:function:: void setValues(const std::vector<dof_no_t> &dofNosLocal, const std::vector<std::array<double,nComponents>> &values, InsertMode petscInsertMode=INSERT_VALUES)
  
  Set values for all components for dofs, after all calls to setValue(s), finishGhostManipulation has to be called to apply the cached changes.
  
  :param const std\:\:vector<dof_no_t> &dofNosLocal: 
  :param const std\:\:vector<std\:\:array<double,nComponents>> &values: 
  :param InsertMode petscInsertMode=INSERT_VALUES: 
  
.. cpp:function:: template<int N>void setValues(const std::array<dof_no_t,N> &dofNosLocal, const std::array<std::array<double,nComponents>,N> &values, InsertMode petscInsertMode=INSERT_VALUES)
  
  Set values for all components for N dofs, after all calls to setValue(s), finishGhostManipulation has to be called to apply the cached changes.
  
  :param const std\:\:array<dof_no_t,N> &dofNosLocal: 
  :param const std\:\:array<std\:\:array<double,nComponents>,N> &values: 
  :param InsertMode petscInsertMode=INSERT_VALUES: 
  
.. cpp:function:: void setValues(int nValues, const std::vector<dof_no_t> &dofNosLocal, const std::vector<std::array<double,nComponents>> &values, InsertMode petscInsertMode=INSERT_VALUES)
  
  Set values for all components for dofs, only nValues values will be set despite potentially more dofNosLocal, after all calls to setValue(s), finishGhostManipulation has to be called to apply the cached changes.
  
  :param int nValues: 
  :param const std\:\:vector<dof_no_t> &dofNosLocal: 
  :param const std\:\:vector<std\:\:array<double,nComponents>> &values: 
  :param InsertMode petscInsertMode=INSERT_VALUES: 
  
.. cpp:function:: void setValue(dof_no_t dofLocalNo, const std::array<double,nComponents> &value, InsertMode petscInsertMode=INSERT_VALUES)
  
  Set a single dof (all components), after all calls to setValue(s), finishGhostManipulation has to be called to apply the cached changes.
  
  :param dof_no_t dofLocalNo: 
  :param const std\:\:array<double,nComponents> &value: 
  :param InsertMode petscInsertMode=INSERT_VALUES: 
  
.. cpp:function:: void setValue(int componentNo, dof_no_t dofLocalNo, double value, InsertMode petscInsertMode)
  
  Set a single dof for a given component, after all calls to setValue(s), finishGhostManipulation has to be called to apply the cached changes.
  
  :param int componentNo: 
  :param dof_no_t dofLocalNo: 
  :param double value: 
  :param InsertMode petscInsertMode: 
  
.. cpp:function:: void setValuesWithGhosts(int componentNo, const std::vector<double> &values, InsertMode petscInsertMode=INSERT_VALUES)
  
  Set values for the specified component for all local dofs, after all calls to setValue(s), finishGhostManipulation has to be called to apply the cached changes.
  
  :param int componentNo: 
  :param const std\:\:vector<double> &values: 
  :param InsertMode petscInsertMode=INSERT_VALUES: 
  
.. cpp:function:: void setValuesWithoutGhosts(int componentNo, const std::vector<double> &values, InsertMode petscInsertMode=INSERT_VALUES)
  
  Set values for the specified component for all local dofs, after all calls to setValue(s), finishGhostManipulation has to be called to apply the cached changes.
  
  :param int componentNo: 
  :param const std\:\:vector<double> &values: 
  :param InsertMode petscInsertMode=INSERT_VALUES: 
  
.. cpp:function:: void setValues(double value)
  
  Set value for all dofs.
  
  :param double value: 
  
.. cpp:function:: void setValuesWithGhosts(const std::vector<std::array<double,nComponents>> &values, InsertMode petscInsertMode=INSERT_VALUES)
  
  Set values for the all component for all local dofs, after all calls to setValue(s), finishGhostManipulation has to be called to apply the cached changes.
  
  :param const std\:\:vector<std\:\:array<double,nComponents>> &values: 
  :param InsertMode petscInsertMode=INSERT_VALUES: 
  
.. cpp:function:: void setValuesWithoutGhosts(const std::vector<std::array<double,nComponents>> &values, InsertMode petscInsertMode=INSERT_VALUES)
  
  Set values for the all component for all local dofs, after all calls to setValue(s), finishGhostManipulation has to be called to apply the cached changes.
  
  :param const std\:\:vector<std\:\:array<double,nComponents>> &values: 
  :param InsertMode petscInsertMode=INSERT_VALUES: 
  
.. cpp:function:: void setValuesWithoutGhosts(const std::array<std::vector<double>,nComponents> &values, InsertMode petscInsertMode=INSERT_VALUES)
  
  Set values for the all component for all local dofs, after all calls to setValue(s), finishGhostManipulation has to be called to apply the cached changes.
  
  :param const std\:\:array<std\:\:vector<double>,nComponents> &values: 
  :param InsertMode petscInsertMode=INSERT_VALUES: 
  
.. cpp:function:: void zeroEntries()
  
  Set value to zero for all dofs.
  
  
.. cpp:function:: void setValues(FieldVariable3D &rhs)
  
  Set values from 3D field variables.
  
  :param FieldVariable3D &rhs: 
  
.. cpp:function:: void setValues(FieldVariable<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>,nComponents> &rhs)
  
  Copy the values from another field variable of the same type.
  
  :param FieldVariable<FunctionSpace\:\:FunctionSpace<Mesh\:\:StructuredDeformableOfDimension<D>,BasisFunctionType>,nComponents> &rhs: 
  
.. cpp:function:: void setValue(dof_no_t dofLocalNo, double value, InsertMode petscInsertMode=INSERT_VALUES)
  
  Only for scalar field variables:  set a single dof (all components) , after all calls to setValue(s), finishGhostManipulation has to be called to apply the cached changes.
  
  :param dof_no_t dofLocalNo: 
  :param double value: 
  :param InsertMode petscInsertMode=INSERT_VALUES: 
  
.. cpp:function:: void setValues(Vec petscVector)
  
  Only for scalar field variables:  set the values from a petsc Vec.
  
  :param Vec petscVector: 
  
.. cpp:function:: void setValues(const std::vector<dof_no_t> &dofNosLocal, std::vector<double> &values, InsertMode petscInsertMode=INSERT_VALUES)
  
  Only for scalar field variables:  set values for the single component for dofs, after all calls to setValue(s), finishGhostManipulation has to be called to apply the cached changes.
  
  :param const std\:\:vector<dof_no_t> &dofNosLocal: 
  :param std\:\:vector<double> &values: 
  :param InsertMode petscInsertMode=INSERT_VALUES: 
  
.. cpp:function:: template<int nValues>void setValues(const std::array<dof_no_t,nValues> dofNosLocal, std::array<double,nValues> values, InsertMode petscInsertMode=INSERT_VALUES)
  
  Only for scalar field variables:  set values for the single component for dofs, after all calls to setValue(s), finishGhostManipulation has to be called to apply the cached changes.
  
  :param const std\:\:array<dof_no_t,nValues> dofNosLocal: 
  :param std\:\:array<double,nValues> values: 
  :param InsertMode petscInsertMode=INSERT_VALUES: 
  
.. cpp:function:: void setValuesWithGhosts(const std::vector<double> &values, InsertMode petscInsertMode=INSERT_VALUES)
  
  Only for scalar field variables:  set values for the single component for all local dofs, after all calls to setValue(s), finishGhostManipulation has to be called to apply the cached changes.
  
  :param const std\:\:vector<double> &values: 
  :param InsertMode petscInsertMode=INSERT_VALUES: 
  

