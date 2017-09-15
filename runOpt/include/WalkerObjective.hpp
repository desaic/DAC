#pragma once
#include "RealFun.hpp"
#include "ElementMesh.hpp"
#include "ElementRegGrid.hpp"
#include "WalkerObjective.hpp"
#include "WorldStair.hpp"

class WalkerObjective : public RealFun{
public:

  WalkerObjective() :
    em(0), world(0), dParam(1e-3),
    nEval(0){}

  ~WalkerObjective(){
    if (em != 0){
      delete em;
    }
  }

  void updateParam(){
  }

  void init(const Eigen::VectorXd & x0);

  ///@brief set simulation to time 0. Restore external forces, x, v, etc.
  void resetSim();

  ///@override
  ///@param x. Two parameters. Width and height.
  void setParam(const Eigen::VectorXd & x0);

  ///@override
  double f();
  Eigen::VectorXd df();

  ///@brief finite element mesh for simulation
  ElementRegGrid * em;

  WorldStair * world;
  WalkerParam walkerParam;

  ConfigFile * conf;
  std::vector<MaterialQuad*> materials;
  float dParam;
  int nEval;
};
