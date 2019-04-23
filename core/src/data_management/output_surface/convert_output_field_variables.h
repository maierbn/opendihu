#pragma once

namespace Data
{

template<typename FieldVariableType>
struct ConvertFieldVariable
{
  typedef FieldVariableType type;

  void convert(const FieldVariableType fieldVariable3D, type &fieldVariable2D, Mesh::face_t face)
  {
    LOG(DEBUG) << "convert field variable: no transformation " << typeid(FieldVariableType).name();
    fieldVariable2D = fieldVariable3D;
  }
};

template<typename BasisFunctionType, int nComponents>
struct ConvertFieldVariable<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>,BasisFunctionType>,nComponents>>>
{
  typedef std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<2>,BasisFunctionType>,nComponents>> type;

  void convert(const std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>,BasisFunctionType>,nComponents>> fieldVariable3D, type &fieldVariable2D, Mesh::face_t face)
  {
    LOG(DEBUG) << "convert field variable: transform shared_ptr<field variable 3D> to shared_ptr<field variable 2D>";
    fieldVariable2D = std::make_shared<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<2>,BasisFunctionType>,nComponents>>>(fieldVariable3D, face);
  }
};

template<typename OutputFieldVariables3D, size_t currentIndex, typename ConvertedTail>
struct ConvertTuple
{
  typedef ConvertTuple<
    OutputFieldVariables3D,        // input tuple of field variables
    currentIndex-1,           // index, starting from which the field variables are already transformed
    decltype(std::tuple_cat(           // tuple of transformed field variables, at the end of the tuple
      std::declval<typename ConvertFieldVariable<std::tuple_element<currentIndex-1,OutputFieldVariables3D>>::type>(),
      std::declval<ConvertedTail>()
    ))
  > ConvertTupleClass;
  typedef typename ConvertTupleClass::type type;

  static void convert(const OutputFieldVariables3D fieldVariables3D, type &fieldVariables2D, Mesh::face_t face)
  {
    LOG(DEBUG) << "convert: currentIndex = " << currentIndex << ", call previous ConvertTuple";

    LOG(DEBUG) << "call convert on field variable";
    ConvertFieldVariable<std::tuple_element<currentIndex,OutputFieldVariables3D>>::convert(std::get<currentIndex>(fieldVariables3D), std::get<currentIndex>(fieldVariables2D), face);

    LOG(DEBUG) << "call convert on previous variables";
    ConvertTupleClass::convert(fieldVariables3D, fieldVariables2D);

  }
};

// recursion end
template<typename OutputFieldVariables3D,typename ConvertedTail>
struct ConvertTuple<OutputFieldVariables3D, 0, ConvertedTail>
{
  typedef ConvertedTail type;

  static void convert(const OutputFieldVariables3D fieldVariables3D, type &fieldVariables2D, Mesh::face_t face)
  {
    LOG(DEBUG) << "convert: recursion end";
    LOG(DEBUG) << "call convert on field variable";
    ConvertFieldVariable<std::tuple_element<0,OutputFieldVariables3D>>::convert(std::get<0>(fieldVariables3D), std::get<0>(fieldVariables2D), face);
  }
};

// interface to be used
template<typename OutputFieldVariables3D>
struct ConvertOutputFieldVariables
{
  typedef ConvertTuple<
    OutputFieldVariables3D,
    std::tuple_size<std::tuple<OutputFieldVariables3D>>::value-1,      // last index of tuple
    std::tuple<    // tuple containing only transformed last element
      typename ConvertFieldVariable   // transformed last element
      <
        std::tuple_element    // last element
        <
          std::tuple_size<std::tuple<OutputFieldVariables3D>>::value-1,
          OutputFieldVariables3D
        >
      >::type
    >
  > ConvertTupleClass;
  typedef typename ConvertTupleClass::type type;   // equals tuple of converted field variables

  static void convert(const OutputFieldVariables3D fieldVariables3D, type &fieldVariables2D, Mesh::face_t face)
  {
    LOG(DEBUG) << "convert: start at interface";
    ConvertTupleClass::convert(fieldVariables3D, fieldVariables2D, face);
  }
};

} // namespace
