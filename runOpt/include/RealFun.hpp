#ifndef REALFUN_HPP
#define REALFUN_HPP

#include <iostream>

#include "Eigen/Dense"

///@brief a real-valued function with n-dimension real input.
///The name ''Objective'' was already used in another file.
///Can be differentiable.
///Can be extended to contain linear or non-linear constraints.

class RealFun{

public:

  Eigen::VectorXd param;

  virtual void init(const Eigen::VectorXd & x0){
    param = x0;
  }

  ///@brief update parameters for computing function value
  ///and gradients. Override to update additional data structures.
  virtual void setParam(const Eigen::VectorXd & x0){
    param = x0;
  }
  ///@brief evaluate function at the currently set parameter.
  virtual double f() = 0;

  ///@brief compute gradient. Not implemented by default.
  ///The return value should have the same dimensions as the parameters.
  virtual Eigen::VectorXd df(){
    std::cout<<"Unimplemented.\n";
    return Eigen::VectorXd(1);
  }

  ///@brief Hessian.
  virtual Eigen::MatrixXd d2f(){
    std::cout<<"Unimplementd.\n";
    return Eigen::MatrixXd(1,1);
  }

  ///@brief update internal parameters.
  virtual void updateParam(){

  }

  virtual ~RealFun(){}
};

#endif 
