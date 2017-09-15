#ifndef NEWTONMIN_HPP
#define NEWTONMIN_HPP

#include "LinSolver.hpp"

#include <Eigen/Sparse>
#include <vector>

class Objective;

///@brief minimize function using Newton's method
///Does not delete any point on destruction.
class NewtonMin{
public:

  NewtonMin():obj(0), maxIters(20), df_rel_tol(1e-3f),
    df_abs_tol(1e-10f), f_tol(1e-6f), residual_tol(1e-11f){}

  ~NewtonMin(){}

  int minimize(std::vector<double> &x);

  Objective * obj;

  int maxIters;

  ///@brief relative tolerance for gradient
  float df_rel_tol;

  ///@brief absolute tolerance for gradient
  float df_abs_tol;

  ///@brief tolerance for change in function value.
  ///Unused in current implementation.
  float f_tol;

  ///@brief tolerance for residual of linear solver.
  float residual_tol;

  LinSolver<double> * solver;
};


#endif
