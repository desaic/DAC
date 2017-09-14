#ifndef QUADRATURE_HPP
#define QUADRATURE_HPP

#include <vector>
#include <Eigen/Dense>

class Quadrature
{
public:
  std::vector<Eigen::Vector3d> x;
  std::vector<float> w;
  Quadrature();
  
  static const Quadrature Gauss2_2D;
  static const Quadrature Gauss4_2D;

  static const Quadrature Gauss2;
  //not a valid quadrature rule.
  static const Quadrature NinePt;
  static const Quadrature Uniform4;
};

#endif
