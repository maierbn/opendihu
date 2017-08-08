#include <iostream>
#include <cstdlib>

#include "opendihu.h"

template<class T1>
class A
{
public:
  T1 t1;
  static void test()
  {
    std::cout<<"A::test()"<<std::endl;
  }
};

class B1
{
public:
  static void test()
  {
    std::cout<<"B1::test()"<<std::endl;
  }
};

class B2
{
public:
  static void test()
  {
    std::cout<<"B2::test()"<<std::endl;
  }
};

template<class T1>
class C
{
public:
  static void test()
  {
    std::cout<<"C::test()"<<std::endl;
  }
};

// partial specialization
template<>
class C<int>
{
public:
  static void test()
  {
    std::cout<<"C<int>::test()"<<std::endl;
  }
};

// different hierarchies possible for partial specialization
template<>
class C<float> : public B1
{
};

template<>
class C<A<int> > : public B2
{
};

int main(int argc, char *argv[])
{
  std::cout<<"template test"<<std::endl;
  A<int>::test();
  C<int>::test();  
  return EXIT_SUCCESS;
}