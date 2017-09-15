#ifndef STEPPERNEWTONDYN_HPP
#define STEPPERNEWTONDYN_HPP

#include "LinSolver.hpp"
#include "Stepper.hpp"

#include <Eigen/Dense>
#include <vector>

class StepperNewtonDyn:public Stepper
{
public:
  StepperNewtonDyn();

  int oneStep();

  float dx_tol;

private:
  ///@param collide. Pass in empty array if no collision
  float compute_dx_sparse(ElementMesh * iMesh,
                          const std::vector<Eigen::Vector3d> & iForces,
                          const std::vector<bool> &collide,
                          std::vector<double> &bb);

  bool m_Init;

};
#endif
