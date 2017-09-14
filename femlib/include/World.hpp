#ifndef WORLDHPP
#define WORLDHPP

#include <vector>
#include <thread>
#include <Eigen/Dense>
#include "Box.hpp"
#include "CCDTrigMesh.hpp"
#include "Contact.hpp"
#include "LinSolver.hpp"

class ElementMesh;
class Stepper;
class TrigMesh;
class Material;
class Quadrature;

class World{

public:
  
  enum SimStage{SIM_LOADING, SIM_LAUNCHING, SIM_FLIGHT};

  World();
  ~World();

  ///@brief deleted by world destructor.
  std::vector<Stepper*> stepperList;

  ///@brief deleted by world destructor.
  std::vector<ElementMesh*> em_;
  std::vector<LinSolver<double>* > solvers;
  std::vector<TrigMesh*> trigm;
  CCDTrigMesh ccd;
  std::vector<Material*>materials;
  const Quadrature * quadrature;
  void launchStepperThread();
  void loop();
  void saveMeshState();
  void resolveCollisionVer0(std::vector<Contact> & collision);
  void resolveCollisionVer1(std::vector<Contact> & collision);

  //configurations

  ///@brief time step to use when in initial contact
  float   contact_h;

  int contactModel;

  bool handle_contact;

  ///@brief time step size for ballistic motion
  float   flight_h;

  ///@brief maximum force allowed before transitioning into rigid body sim
  float   flight_fmax;

  Box box;

  ///@brief 0 for nothing. 1 for floor. 2 for box.
  int collideWith;

  Eigen::Vector3d boxContactNormal;

  ///@brief stop sim if there is collision.
  bool   terminateOnCollision;

  //////////
  //time varying simulation states.

  ///@brief mesh being processed
  int   activeMesh;

  ///@brief number of time steps taken.
  int     frameCnt;

  ///@brief highest point of jumper.
  double height;

  ///@brief virtual world time spent in loading in seconds
  float loading_t;

  ///@brief is it during force loading
  SimStage stage;

  ///@brief pointer to active stepper.
  Stepper * stepper;

  ///@brief which time stepper is in use.
  int stepperIdx;

  int timeStep;

  ///@brief total rotation.
  double totalAng;
  
  ///@brief directory for simulation state position and velocity
  std::string   simDirName;

  ///@brief save every n frames.
  int   frameSkip;

  ///@brief save position and velocity to files?
  bool   savePos;

};
#endif
