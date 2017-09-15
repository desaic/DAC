#ifndef STEPPERNEWMARK_HPP
#define STEPPERNEWMARK_HPP

#include "Stepper.hpp"
#include "LinSolver.hpp"
#include "NewtonMin.hpp"
#include <Eigen/Dense>
#include <vector>
#include <fstream>

class ElementMesh;

class StepperNewmark: public Stepper
{
public:
  StepperNewmark();
  ~StepperNewmark();
   void init();

    int oneStep();

  float getEnergy();

  int runNewton(std::vector<double> &x);

  std::vector<double> getGrad(bool constrained=true);
  std::vector<double> getRHS();
  Eigen::SparseMatrix<double> getHessian();

  void setx(double *x);

  void resolveCollision(std::vector<Contact> & collision);

  //void resolveCollision2D(std::vector<bool> & collision);

  int resolveCollision3D_old(std::vector<Contact> & collision);
  void resolveCollision3D(std::vector<Contact> & collision);

  void updateStiffness();

  void updateMb();

//private:
  bool update_stiffness;
  Eigen::VectorXd delta;//, delta0;
  std::vector<Eigen::Vector3d> x0,v0,Mb;

  ///@brief damping matrix
  Eigen::SparseMatrix<double> DampMat,K;

  bool m_Init;
  bool useOldContact;
  NewtonMin * newton;
};

#endif
