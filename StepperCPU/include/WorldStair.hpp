#ifndef WORLD_STAIR_HPP
#define WORLD_STAIR_HPP

#include "ConfigFile.hpp"
#include "ElementRegGrid.hpp"
#include "LinSolver.hpp"
#include "MaterialQuad.hpp"
#include "Stepper.hpp"
#include "TrigMesh.hpp"
#include <vector>
#include <Eigen/Dense>

struct Contact;

class WalkerParam
{
public:
  
  //x y z size.
  float legSize[3];
  
  //size of middle portion of body
  float bodySize[3];
  
  //target element size. May be off by a little.
  float dx;
  
  float minBodySize[3];

  WalkerParam() :dx(0.01f){
    legSize[0] = dx;
    legSize[1] = 3*dx;
    legSize[2] = dx;

    bodySize[0] = 3* dx;
    bodySize[1] = dx;
    bodySize[2] = 3 * dx;

    minBodySize[0] = 2 * dx;
    minBodySize[1] = 0.1 * dx;
    minBodySize[2] = 2 * dx;

  }
};

//simulation for catapults
class WorldStair
{
public :

  WorldStair();
  ~WorldStair();

  ///@brief deleted by world destructor.
  std::vector<ElementMesh*> em_;
  std::vector<LinSolver<double>* > solvers;
  std::vector<Material*>materials;
  const Quadrature * quadrature;
  void saveMeshState();

  int  resolveCollision(std::vector<Contact> & contact);
  void loop();
  bool finishSim();
  //@brief initialize linear solver data structures.
  //computes numerical factorization.
  void initSolvers();
  ///@brief set everything to starting configuration.
  void resetSim();
  
  ///@brief time step
  float   dt;
  float   forceTime;
  float dtheta;
  std::vector<std::vector<Eigen::Vector3d> > fe;
  bool hasInitPos;
  float stairHeight;
  float stairHeightStart;
  float stairWidth;
  float wallDist;

  /////////////////////////////////
  //time varying simulation states.
  /////////////////////////////////
  
  ///@brief mesh being processed
  int   activeMesh;

  std::vector<Contact> contact;

  ///@brief number of time steps taken.
  int     frameCnt;

  //center of mass at beginning of sim.
  Eigen::Vector3d centerStart;
  Eigen::Vector3d centerEnd;
  // time step at which the sim finished.
  float finish_time;
  bool past_finish_line;
  //speed is calculated as displacement in x y direction divide by finish_time.
  float overall_speed;

  ///@brief pointer to active stepper.
  std::vector<Stepper *> stepperList;
  Stepper * stepper;
  //=========================
  //visualization and logging
  //=========================
  ///@brief directory for simulation state position and velocity
  std::string   simDirName;

  ///@brief save every n frames.
  int   frameSkip;

  ///@brief save position and velocity to files?
  bool   savePos;
};

int loadStairScene(WorldStair * world, const ConfigFile & conf);

void makeWalker(ElementRegGrid * em, WalkerParam * param);

#endif
