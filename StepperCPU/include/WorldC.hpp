#ifndef WORLD_C_HPP
#define WORLD_C_HPP

#include "Box.hpp"
#include "ConfigFile.hpp"
#include "LinSolver.hpp"
#include "MaterialQuad.hpp"
#include "Stepper.hpp"
#include "TrigMesh.hpp"
#include <vector>
#include <Eigen/Dense>

struct Contact;
class StepperStatic;

struct PositionOnElement
{
  ///@brief element index.
  int eidx;
  ///@brief bilinear weights
  float bw[4];
  PositionOnElement() :eidx(-1){}
};

//simulation for catapults
class WorldC
{
public :
  WorldC();
  ~WorldC();

  ///@brief deleted by world destructor.
  std::vector<Stepper*> stepperList;
  StepperStatic * staticSim;
  ///@brief deleted by world destructor.
  std::vector<ElementMesh*> em_;
  std::vector<LinSolver<double>* > solvers;
  std::vector<TrigMesh*> trigm;
  std::vector<Material*>materials;
  const Quadrature * quadrature;
  void saveMeshState();

  void resolveCollision(std::vector<Contact> & contact);
  void loop();
  void solveStatic();

  //@brief initialize linear solver data structures.
  //computes numerical factorization.
  void initSolvers();
  ///@brief set everything to starting configuration.
  void resetSim();

  //////////////////////////////////////////////////////
  ///sim configurations that does not change over time.
  /////////////////////////////////////////////////////

  Box catapultBox;
  
  ///@brief where to place the cube relative to the
  ///left, top, and back side of the catapult.
  Eigen::Vector3d cubeLocation;
  
  ///@brief time step to use when in initial contact
  float   contact_h;

  int contactModel;

  ///@brief initial external forces.
  std::vector<std::vector<Eigen::Vector3d> > fe;

  bool handle_contact;

  bool hasInitPos;

  ///@brief time step size for ballistic motion
  float   flight_h;
  
  float launchAreaSize;

  float targetHeight;

  ///@brief stop sim if there is collision.
  bool   terminateOnCollision;

  ///@brief distance of wall the projectile hits.
  float wallDist;

  /////////////////////////////////
  //time varying simulation states.
  /////////////////////////////////
  
  ///@brief mesh being processed
  int   activeMesh;

  PositionOnElement cubeFixture;

  std::vector<Contact> contact;

  ///@brief number of time steps taken.
  int     frameCnt;
  
  typedef enum{LAUNCH, FLIGHT}SimStage;

  SimStage stage;

  ///@brief pointer to active stepper.
  Stepper * stepper;

  ///@brief which time stepper is in use.
  int stepperIdx;

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
int loadMaterials(WorldC * world, const ConfigFile & conf);
int loadCScene(WorldC * world, const ConfigFile & conf);

#endif
