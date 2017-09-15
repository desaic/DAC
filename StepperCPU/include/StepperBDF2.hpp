#ifndef STEPPERBDF2_HPP
#define STEPPERBDF2_HPP

#include "Stepper.hpp"
#include "LinSolver.hpp"
#include "NewtonMin.hpp"
#include <Eigen/Dense>
#include <vector>
#include <fstream>

class ElementMesh;

class StepperBDF2: public Stepper
{
public:
  StepperBDF2();
  ~StepperBDF2();
   void init();

    int oneStep();

  float getEnergy();

   int runNewton(std::vector<double> &x);

  std::vector<double> getGrad(bool constrained=true);

  Eigen::SparseMatrix<double> getHessian();

  void setx(double *x);

  void resolveCollision(std::vector<bool> & collision);

  void resolveCollision2D(std::vector<bool> & collision);

  void resolveCollision3D(std::vector<bool> & collision);

  void updateStiffness();

  //update explicit terms
  void update_p0();

private:
  bool update_stiffness;
  Eigen::VectorXd delta;
  std::vector<Eigen::Vector3d> x0,x1,v0,v1;
  Eigen::VectorXd p0;

  ///@brief damping matrix
  Eigen::SparseMatrix<double> D,K;

  bool m_Init;

  NewtonMin * newton;
};

#endif
