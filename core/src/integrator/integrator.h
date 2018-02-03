#pragma once

namespace Integrator
{

class Integrator
{
public:
};

/**
 * Integrator class that can be used if no integrator is necessary, because 
 * the integral is solve analytically. This is the case if Mesh::RegularFixed is used.
 * Those meshes use stencils.
 */
class None
{
public:
  static constexpr int numberEvaluations(){return 0;}
};

};
