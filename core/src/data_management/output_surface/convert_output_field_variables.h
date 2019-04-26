#pragma once

namespace Data
{

template<typename FieldVariableType>
struct ConvertFieldVariable
{
  typedef FieldVariableType type;

  static void convert(const FieldVariableType fieldVariable3D, type &fieldVariable2D, Mesh::face_t face, bool &ownRankInvolvedInOutput)
  {
    LOG(DEBUG) << "convert field variable \"" << fieldVariable3D->name() << "\": no transformation " << typeid(FieldVariableType).name();
    fieldVariable2D = fieldVariable3D;
  }
};

template<typename BasisFunctionType, int nComponents>
struct ConvertFieldVariable<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>,BasisFunctionType>,nComponents>>>
{
  typedef std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<2>,BasisFunctionType>,nComponents>> type;

  static void convert(const std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>,BasisFunctionType>,nComponents>> fieldVariable3D, type &fieldVariable2D, Mesh::face_t face, bool &ownRankInvolvedInOutput)
  {
    LOG(DEBUG) << "convert field variable \"" << fieldVariable3D->name() << "\": transform shared_ptr<field variable 3D> to shared_ptr<field variable 2D>";
    if (!fieldVariable2D)
    {
      fieldVariable2D = std::make_shared<FieldVariable::FieldVariable<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<2>,BasisFunctionType>,nComponents>>(*fieldVariable3D, face, ownRankInvolvedInOutput);
    }
    else
    {
      fieldVariable2D->setValues(*fieldVariable3D);
    }
  }
};

template<typename OutputFieldVariables3D, size_t currentIndex, typename ConvertedTail>
struct ConvertTuple
{
  typedef ConvertTuple<
    OutputFieldVariables3D,        // input tuple of field variables
    currentIndex-1,           // index, starting from which the field variables are already transformed
    decltype(std::tuple_cat(           // tuple of transformed field variables, at the end of the tuple
      std::declval<std::tuple<
          typename ConvertFieldVariable<typename std::tuple_element<currentIndex-1,OutputFieldVariables3D>::type>::type
          >>(),
      std::declval<ConvertedTail>()
    ))
  > ConvertTupleClass;
  typedef typename ConvertTupleClass::type type;

  // convert a tuple of field variables from 3D function spaces to 2D function spaces
  static void convert(const OutputFieldVariables3D fieldVariables3D, type &fieldVariables2D, Mesh::face_t face, bool &ownRankInvolvedInOutput)
  {
    LOG(DEBUG) << "convert: currentIndex = " << currentIndex << ", call previous ConvertTuple";

    LOG(DEBUG) << "call convert on field variable";
    ConvertFieldVariable<typename std::tuple_element<currentIndex,OutputFieldVariables3D>::type>::convert(
      std::get<currentIndex>(fieldVariables3D), std::get<currentIndex>(fieldVariables2D), face, ownRankInvolvedInOutput
    );

    LOG(DEBUG) << "call convert on previous variables";
    ConvertTupleClass::convert(fieldVariables3D, fieldVariables2D, face, ownRankInvolvedInOutput);

  }
};

// recursion end
template<typename OutputFieldVariables3D,typename ConvertedTail>
struct ConvertTuple<OutputFieldVariables3D, 0, ConvertedTail>
{
  typedef ConvertedTail type;

  static void convert(const OutputFieldVariables3D fieldVariables3D, type &fieldVariables2D, Mesh::face_t face, bool &ownRankInvolvedInOutput)
  {
    LOG(DEBUG) << "convert: recursion end";
    LOG(DEBUG) << "call convert on field variable";
    ConvertFieldVariable<typename std::tuple_element<0,OutputFieldVariables3D>::type>::convert(
      std::get<0>(fieldVariables3D), std::get<0>(fieldVariables2D), face, ownRankInvolvedInOutput
    );
  }
};

// interface to be used
template<typename OutputFieldVariables3D>
struct ConvertOutputFieldVariables
{
  typedef ConvertTuple<
    OutputFieldVariables3D,
    std::tuple_size<OutputFieldVariables3D>::value-1,      // last index of tuple
    std::tuple<    // tuple containing only transformed last element
      typename ConvertFieldVariable   // transformed last element
      <
        typename std::tuple_element    // last element
        <
          std::tuple_size<OutputFieldVariables3D>::value-1,
          OutputFieldVariables3D
        >::type
      >::type
    >
  > ConvertTupleClass;
  typedef typename ConvertTupleClass::type type;   // equals tuple of converted field variables

  typedef typename std::tuple_element<0,type>::type::element_type::FunctionSpace FunctionSpaceFirstFieldVariable;  // the function space type of the first field variable

  // convert a tuple of field variables from 3D function spaces to 2D function spaces by taking the surface at the given face
  static void convert(const OutputFieldVariables3D fieldVariables3D, type &fieldVariables2D, Mesh::face_t face, bool &ownRankInvolvedInOutput)
  {
    LOG(DEBUG) << "convert: start 3D: " << typeid(OutputFieldVariables3D).name() << ", n field variables: " << std::tuple_size<std::tuple<OutputFieldVariables3D>>::value;
    ConvertTupleClass::convert(fieldVariables3D, fieldVariables2D, face, ownRankInvolvedInOutput);
  }

  // extract the function space of the first field variable
  static std::shared_ptr<FunctionSpaceFirstFieldVariable> getFunctionSpaceFirstFieldVariable(type &fieldVariables2D)
  {
    return std::get<0>(fieldVariables2D)->functionSpace();
  }
};

} // namespace
