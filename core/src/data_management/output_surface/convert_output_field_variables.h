#pragma once

namespace Data
{

template<typename FieldVariableType>
struct ConvertFieldVariable
{
  typedef FieldVariableType type;

  static void convert(const FieldVariableType fieldVariable3D, type &fieldVariable2D, std::vector<Mesh::face_t> &faces, bool &ownRankInvolvedInOutput)
  {
    VLOG(1) << "convert field variable \"" << fieldVariable3D->name() << "\": no transformation " << StringUtility::demangle(typeid(FieldVariableType).name());
    fieldVariable2D = fieldVariable3D;
  }
};

template<typename BasisFunctionType, int nComponents>
struct ConvertFieldVariable<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>,BasisFunctionType>,nComponents>>>
{
  typedef std::vector<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<2>,BasisFunctionType>,nComponents>>> type;

  static void convert(const std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>,BasisFunctionType>,nComponents>> fieldVariable3D,
                      type &fieldVariables2D, std::vector<Mesh::face_t> &faces, bool &ownRankInvolvedInOutput)
  {
    using FieldVariable2D = FieldVariable::FieldVariable<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<2>,BasisFunctionType>,nComponents>;

    VLOG(1) << "convert field variable \"" << fieldVariable3D->name() << "\": transform shared_ptr<field variable 3D> to vector<shared_ptr<field variable 2D>>";
    if (fieldVariables2D.empty())
    {
      // fill vector with field variables
      for (Mesh::face_t face : faces)
      {
        std::shared_ptr<FieldVariable2D> fieldVariable2D = std::make_shared<FieldVariable2D>(*fieldVariable3D, face, ownRankInvolvedInOutput);
        fieldVariables2D.push_back(fieldVariable2D);
      }
    }
    else
    {
      assert(fieldVariables2D.size() == faces.size());

      // fill vector with field variables
      for (int i = 0; i < faces.size(); i++)
      {
        // this just updates the values, which dofs to use is already stored inside the field variable
        fieldVariables2D[i]->setValues(*fieldVariable3D);
      }
    }
  }
};

template<typename BasisFunctionType, int nComponents>
struct ConvertFieldVariable<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<3>,BasisFunctionType>,nComponents>>>
{
  typedef std::vector<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<2>,BasisFunctionType>,nComponents>>> type;

  static void convert(const std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<3>,BasisFunctionType>,nComponents>> fieldVariable3D,
                      type &fieldVariables2D, std::vector<Mesh::face_t> &faces, bool &ownRankInvolvedInOutput)
  {
    using FieldVariable2D = FieldVariable::FieldVariable<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<2>,BasisFunctionType>,nComponents>;

    VLOG(1) << "convert composite field variable \"" << fieldVariable3D->name() << "\": transform shared_ptr<field variable 3D> to vector<shared_ptr<field variable 2D>>";
    if (fieldVariables2D.empty())
    {
      ownRankInvolvedInOutput = false;
      bool ownRankInvolvedInOutputForFace = false;

      // fill vector with field variables
      int i = 0;
      for (Mesh::face_t face : faces)
      {
        std::shared_ptr<FieldVariable2D> fieldVariable2D = std::make_shared<FieldVariable2D>(*fieldVariable3D, face, ownRankInvolvedInOutputForFace);
        fieldVariables2D.push_back(fieldVariable2D);

        if (ownRankInvolvedInOutputForFace)
          ownRankInvolvedInOutput = true;

        LOG(DEBUG) <<  "i=" << i << "/" << faces.size() << ", initialize field variable \"" << fieldVariable2D->name() << "\", values: " << fieldVariable2D->partitionedPetscVec();
        i++;
      }
    }
    else
    {
      assert(fieldVariables2D.size() == faces.size());

      // fill vector with field variables
      for (int i = 0; i < faces.size(); i++)
      {
        LOG(DEBUG) << "i=" << i << "/" << faces.size() << ", use field variable \""
          << fieldVariables2D[i]->name() << "\", values: " << fieldVariables2D[i]->partitionedPetscVec();

        // this just updates the values, which dofs to use is already stored inside the field variable
        fieldVariables2D[i]->setValues(*fieldVariable3D);
      }
    }
  }
};

template<typename FieldVariableType>
struct ConvertFieldVariable<std::vector<FieldVariableType>>
{
  typedef std::vector<typename ConvertFieldVariable<FieldVariableType>::type> type;

  static void convert(const std::vector<FieldVariableType> &fieldVariable3D,
                      type &fieldVariable2D, std::vector<Mesh::face_t> &faces, bool &ownRankInvolvedInOutput)
  {
    fieldVariable2D.resize(fieldVariable3D.size());

    for (int i = 0; i < fieldVariable3D.size(); i++)
    {
      ConvertFieldVariable<FieldVariableType>::convert(fieldVariable3D[i], fieldVariable2D[i], faces, ownRankInvolvedInOutput);
    }
  }
};

template<typename FieldVariablesForOutputWriter3D, size_t currentIndex, typename ConvertedTail>
struct ConvertTuple
{
  typedef ConvertTuple<
    FieldVariablesForOutputWriter3D,        // input tuple of field variables
    currentIndex-1,           // index, starting from which the field variables are already transformed
    decltype(std::tuple_cat(           // tuple of transformed field variables, at the end of the tuple
      std::declval<std::tuple<
          typename ConvertFieldVariable<typename std::tuple_element<currentIndex-1,FieldVariablesForOutputWriter3D>::type>::type
          >>(),
      std::declval<ConvertedTail>()
    ))
  > ConvertTupleClass;
  typedef typename ConvertTupleClass::type type;

  // convert a tuple of field variables from 3D function spaces to 2D function spaces
  static void convert(const FieldVariablesForOutputWriter3D fieldVariables3D,
                      type &fieldVariables2D, std::vector<Mesh::face_t> &faces, bool &ownRankInvolvedInOutput)
  {
    VLOG(1) << "convert: currentIndex = " << currentIndex << ", call previous ConvertTuple, convertedTail: " << StringUtility::demangle(typeid(ConvertedTail).name());

    VLOG(1) << "call convert on field variable " << StringUtility::demangle(typeid(typename std::tuple_element<currentIndex,FieldVariablesForOutputWriter3D>::type).name());
    ConvertFieldVariable<typename std::tuple_element<currentIndex,FieldVariablesForOutputWriter3D>::type>::convert(
      std::get<currentIndex>(fieldVariables3D), std::get<currentIndex>(fieldVariables2D), faces, ownRankInvolvedInOutput
    );

    VLOG(1) << "call convert on previous variables";
    ConvertTupleClass::convert(fieldVariables3D, fieldVariables2D, faces, ownRankInvolvedInOutput);

  }
};

// recursion end
template<typename FieldVariablesForOutputWriter3D,typename ConvertedTail>
struct ConvertTuple<FieldVariablesForOutputWriter3D, 0, ConvertedTail>
{
  typedef ConvertedTail type;

  static void convert(const FieldVariablesForOutputWriter3D fieldVariables3D,
                      type &fieldVariables2D, std::vector<Mesh::face_t> &faces, bool &ownRankInvolvedInOutput)
  {
    VLOG(1) << "convert: recursion end, convertedTail: " << StringUtility::demangle(typeid(ConvertedTail).name());
    VLOG(1) << "call convert on field variable " << StringUtility::demangle(typeid(typename std::tuple_element<0,FieldVariablesForOutputWriter3D>::type).name());
    ConvertFieldVariable<typename std::tuple_element<0,FieldVariablesForOutputWriter3D>::type>::convert(
      std::get<0>(fieldVariables3D), std::get<0>(fieldVariables2D), faces, ownRankInvolvedInOutput
    );
  }
};

// interface to be used
template<typename FieldVariablesForOutputWriter3D>
struct ConvertFieldVariablesForOutputWriter
{
  typedef ConvertTuple<
    FieldVariablesForOutputWriter3D,
    std::tuple_size<FieldVariablesForOutputWriter3D>::value-1,      // last index of tuple
    std::tuple<    // tuple containing only transformed last element
      typename ConvertFieldVariable   // transformed last element
      <
        typename std::tuple_element    // last element
        <
          std::tuple_size<FieldVariablesForOutputWriter3D>::value-1,
          FieldVariablesForOutputWriter3D
        >::type
      >::type
    >
  > ConvertTupleClass;
  typedef typename ConvertTupleClass::type type;   // equals tuple of converted field variables

  typedef typename std::tuple_element<0,type>::type::value_type::element_type FirstFieldVariable;  // the first field variable
  typedef typename FirstFieldVariable::FunctionSpace FunctionSpaceFirstFieldVariable;  // the function space type of the first field variable
  typedef typename std::tuple_element<2,type>::type::value_type::element_type SecondFieldVariable;  // the first field variable

  // convert a tuple of field variables from 3D function spaces to 2D function spaces by taking the surface at the given face
  static void convert(const FieldVariablesForOutputWriter3D fieldVariables3D,
                      type &fieldVariables2D, std::vector<Mesh::face_t> &faces, bool &ownRankInvolvedInOutput)
  {
    VLOG(1) << "convert: start 3D: " << StringUtility::demangle(typeid(FieldVariablesForOutputWriter3D).name())
      << ", n field variables: " << std::tuple_size<std::tuple<FieldVariablesForOutputWriter3D>>::value;
    ConvertTupleClass::convert(fieldVariables3D, fieldVariables2D, faces, ownRankInvolvedInOutput);
  }

  // extract the function space of the first field variable
  static std::shared_ptr<FunctionSpaceFirstFieldVariable> getFunctionSpaceFirstFieldVariable(type &fieldVariables2D)
  {
    VLOG(1) << "getFunctionSpaceFirstFieldVariable, type of fieldVariable2D: " << StringUtility::demangle(typeid(type).name());
    int n2DFunctionSpaces = std::get<0>(fieldVariables2D).size();
    if (n2DFunctionSpaces > 0)
      return std::get<0>(fieldVariables2D)[0]->functionSpace();

    return nullptr;
  }

  // extract all the function spaces for the different faces of the first field variable
  static void getFunctionSpacesFirstFieldVariable(type &fieldVariables2D, std::vector<std::shared_ptr<FunctionSpaceFirstFieldVariable>> &functionSpaces)
  {
    int n2DFunctionSpaces = std::get<0>(fieldVariables2D).size();
    functionSpaces.resize(n2DFunctionSpaces);

    for (int i = 0; i < n2DFunctionSpaces; i++)
    {
      functionSpaces[i] = std::get<0>(fieldVariables2D)[i]->functionSpace();
    }
  }
};

} // namespace
