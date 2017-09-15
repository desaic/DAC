#ifndef STEPPERLINNEWMARK_HPP
#define STEPPERLINNEWMARK_HPP

#include "Stepper.hpp"
#include "LinSolver.hpp"
#include <Eigen/Dense>
#include <vector>

class ElementMesh;

class StepperLinNewmark: public Stepper
{
public:
  StepperLinNewmark();

    int oneStep();

private:
  ///@param collide. Pass in empty array if no collision
  double solve_acceleration(const std::vector<Eigen::Vector3d> & forces,
                          const std::vector<bool> & collide,
                          std::vector<double> & bb);

  ///@brief acceleration for next time step.
  std::vector<double>  anext;

  bool m_Init;
};

#endif
