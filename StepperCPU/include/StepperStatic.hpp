#ifndef STEPPER_STATIC_HPP
#define STEPPER_STATIC_HPP

#include "Stepper.hpp"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>

class NewtonMin;

class StepperStatic:public Stepper
{
public:
  StepperStatic();
  ~StepperStatic();
  void init();
  int oneStep();

  float getEnergy();

  std::vector<double> getGrad(bool constrained=true);

  Eigen::SparseMatrix<double> getHessian();

  void setx(double * x);

  //convergence criteria
  double force_relative_tol;
  double force_abs_tol;
  double force_norm0;

  //intermediate sim states
  std::vector<int> fixed0;

  //equality constraints.
  //each vector contains a list of dofs that needs to be equal to each other.
  std::vector<std::vector<int> > CEq;

private:
  std::vector<double> delta;
  std::vector<Eigen::Vector3d> x0;
  bool m_Init;
  NewtonMin * newton;
};
#endif
