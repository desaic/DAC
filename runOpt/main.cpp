#include <iostream>
#include <fstream>

#include "ConfigUtil.hpp"
#include "FileUtil.hpp"
#include "EigenUtil.hpp"
#include "JumperObjective.hpp"
#include "runLBFGS.hpp"
#include "Material.hpp"
#include "MaterialQuad.hpp"
#include "StrainEne.hpp"
#include "StepperNewmark.hpp"
#include "CatapultObjective.hpp"
#include "WorldC.hpp"
#include "IpIpoptApplication.hpp"
#include "IpSolveStatistics.hpp"
#include "IpInterface.hpp"
#include "TestIpopt.hpp"
#include "WalkerObjective.hpp"

int runIpopt(const ConfigFile & config, RealFun & objective)
{
  ///@TODO: bounds for jumper opt.
  //float lb = -0.007;
  //float ub = 0.007;

  IpInterface * j = new IpInterface();
  j->logger = new std::ofstream("optlog.txt");
  j->x_L.resize(objective.param.size(), 0);

  //bounds for walker.
  //body size.
  j->x_L[0] = objective.param[0];
  j->x_L[1] = 0.2 * objective.param[1];
  j->x_L[2] = objective.param[2];
  //leg size.
  j->x_L[3] = 0.2 * objective.param[3];
  j->x_L[4] = 0.5 * objective.param[4];
  j->x_L[5] = 0.2 * objective.param[5];
  
  j->x_U.resize(objective.param.size(), 0);

  j->x_U[0] = 2 * objective.param[0];
  j->x_U[1] = 2 * objective.param[1];
  j->x_U[2] = 2 * objective.param[2];
  //leg size.
  j->x_U[3] = 2 * objective.param[3];
  j->x_U[4] = 2 * objective.param[4];
  j->x_U[5] = 2 * objective.param[5];

  j->objective = &objective;
  j->init();
  j->run();

  return 0;
}

int main(int argc, char * argv[])
{
  if(argc<2){
    testIpInterface();
    std::cout<<"Opt config_file \n"<<"\n";
    return 0;
  }
  ConfigFile * config = new ConfigFile();
  int status = config->load(argv[1]);
  if(status<0){
    return 0;
  }

  std::string scene = config->getString("scene");
  if (scene == "jumper"){
    JumperObjective objective;
    World * world = new World();
    loadScene(world, *config);
    objective.world = world;
    objective.conf = config;
    int N = 2;
    Eigen::VectorXd param = Eigen::VectorXd::Zero(N);
    for (int ii = 0; ii<param.rows(); ii++){
      param[ii] = 0;
    }
    //restore previous progress
    //param[0] = -0.00477471;
    //param[1] = 0.00636058;
    objective.init(param);
    //  objective.setParam(param);
    runIpopt(*config, objective);
  }
  else if (scene == "stair"){
    WalkerObjective objective;
    WorldStair * world = new WorldStair();
    loadStairScene(world, *config);
    int N = 6;
    Eigen::VectorXd param = Eigen::VectorXd::Zero(N);
    param[0] = 0.02;
    param[1] = 0.01;
    param[2] = 0.02;
    param[3] = 0.005;
    param[4] = 0.02;
    param[5] = 0.01;
    objective.em = (ElementRegGrid*)world->em_[0];
    objective.world = world;
    objective.init(param);
    objective.param = param;
    runIpopt(*config, objective);
  }
  else{
    CatapultObjective objective;
    config->getFloat("dparam", objective.dParam);
    WorldC * world = new WorldC();

    loadCScene(world, *config);
    //optimization does not save meshes.
    world->savePos = false;
    objective.world = world;
    objective.conf = config;
    int N = 2;
    Eigen::VectorXd param = Eigen::VectorXd::Zero(N);
    for (int ii = 0; ii<param.rows(); ii++){
      param[ii] = 0;
    }
    objective.init(param);
    //  objective.setParam(param);
    //  int ret = runLBFGS(config, objective);
    runIpopt(*config, objective);
  }

  ////debug
//  param[0] = 0;
//  param[1] = 0;
//  objective.setParam(param);
//  objective.resetSim();
//  objective.stepper->nSteps = 1500;
//  objective.stepper->step();
//  objective.stepper->nSteps = 10000;
//  param[0] =  0.001261867;
//  param[1] = -0.000455077;
//  objective.setParam(param);

//  for(int ii = 0; ii<param.rows(); ii++){
//    param[ii] = 1e-3;
//  }
//  for(int s = 0; s<10000; s++){
//  objective.setParam(param);
//  objective.f();
//  }
  return 0;
}
