#pragma once

#include "discretizable_in_time/discretizable_in_time.h"

namespace ModelOrderReduction
{

/** A class to distinguish different operations of the POD operation, not yet used
  */
struct LinearPart{};

/** A class that implements the POD model order reduction technique.
 */
template<typename DiscretizableInTimeType>
class PODBase : public DiscretizableInTime
{
public:
  //! constructor
  PODBase(DihuContext context);

  //! initialize timestepping
  void initialize();

  //! timestepping rhs function f of equation u_t = f(u,t)
  virtual void evaluateTimesteppingRightHandSide(Vec &input, Vec &output, int timeStepNo, double currentTime);

  //! get the number of degrees of freedom per node which is 1 by default
  int nComponentsNode();

  //! set initial values and return true or don't do anything and return false
  bool setInitialValues(Vec &initialValues);

  //! return whether the object has a specified mesh type and is not independent of the mesh type
  bool knowsMeshType();

  //! return the used mesh
  std::shared_ptr<Mesh::Mesh> mesh();

protected:

  const DihuContext context_;    ///< object that contains the python config for the current context and the global singletons meshManager and solverManager
  DiscretizableInTimeType problem_;   ///< the DiscretizableInTime object that is managed by this class
};

/** The POD class that is used from main function
 */
template<typename DiscretizableInTimeType, typename PartType>
class POD{};

/** Specialization for "LinearPart"
 */
template<typename DiscretizableInTimeType>
class POD<DiscretizableInTimeType, LinearPart>:
  public PODBase<DiscretizableInTimeType>
{
public:
  typedef typename DiscretizableInTimeType::FunctionSpace FunctionSpace;

  //! use constructor of base class
  using PODBase<DiscretizableInTimeType>::PODBase;

  //! timestepping rhs function f of equation u_t = f(u,t)
  void evaluateTimesteppingRightHandSide(Vec &input, Vec &output, int timeStepNo, double currentTime);
};

};  // namespace

#include "model_order_reduction/pod.tpp"