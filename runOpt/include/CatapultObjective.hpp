#ifndef CATAPULT_OBJECTIVE_HPP
#define CATAPULT_OBJECTIVE_HPP
#include "ConfigFile.hpp"
#include "RealFun.hpp"
#include "StretchMesh.hpp"
#include <Eigen/Dense>
#include <vector>
class WorldC;

class CatapultObjective: public RealFun{
    
public:
  CatapultObjective();

  void init(const Eigen::VectorXd & x0);

  ///@override
  ///@param x. Two parameters. Width and height.
  void setParam(const Eigen::VectorXd & x0);

  ///@override
  double f();
  Eigen::VectorXd df();

  ///@brief step size for numerical differencing.
  float dParam;
  WorldC * world;
  ConfigFile * conf;
  StretchOpts opts;
  
  ///@brief unstretched catapult rest pose.
  std::vector<Eigen::Vector3d> X0;
};

#endif
