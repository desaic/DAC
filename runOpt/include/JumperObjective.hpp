#ifndef JUMPEROBJECTIVE
#define JUMPEROBJECTIVE

#include "ConfigFile.hpp"
#include "ElementMesh.hpp"
#include "RealFun.hpp"
#include "StretchMesh.hpp"

#include <vector>

class MaterialQuad;
class Mesh;
class StepperNewmark;
class World;

///@brief states used to restore sim.
struct WorldState{
  std::vector<Eigen::Vector3d> fe;

  std::vector<int> fixedDof;

  std::vector<Eigen::Vector3d> lb;

  std::vector<Eigen::Vector3d> X;
};

class JumperObjective: public RealFun{

public:

  JumperObjective() :
    em(0), world(0), dParam(2e-4),
    nEval(0){}

  ~JumperObjective(){
    if(em!=0){
      delete em;
    }
  }
  
  std::vector<double> lambda;
  std::vector<double> rho;
  std::vector<double> g;

  void init(const Eigen::VectorXd & x0);
  
  void updateParam(){
    std::cout << "update lambda: \n";
    for (int i = 0; i<lambda.size(); i++){
      lambda[i] = rho[i] * std::max(g[i] + lambda[i] / rho[i], 0.0);
      std::cout << lambda[i] << " ";
    }
    std::cout << "\n";
  }

  ///@brief set simulation to time 0. Restore external forces, x, v, etc.
  void resetSim();

  ///@override
  ///@param x. Two parameters. Width and height.
  void setParam( const Eigen::VectorXd & x0);

  ///@override
  double f();
  Eigen::VectorXd df();
  ///@brief finite element mesh for simulation
  ElementMesh * em;

  World * world;

  ConfigFile * conf;

  WorldState state;

  StretchOpts opts;

  bool updateG;
  //central difference step size.
  double dParam;
  ///@brief number of calls made to f().
  int nEval;

  std::vector<MaterialQuad*> materials;
};

#endif
