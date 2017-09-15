#ifndef WORLD_COMPRESS_HPP
#define WORLD_COMPRESS_HPP

#include "ConfigFile.hpp"
#include "LinSolver.hpp"
#include "Quadrature.hpp"

class ElementMesh;
class StepperStatic;

class WorldCompress{
public:
  WorldCompress();
  ~WorldCompress();
  
  ElementMesh* em_;
  std::vector<Material*>materials;
  LinSolver<double>* solver;
  const Quadrature * quadrature;
 
  //@brief initialize linear solver data structures.
  //computes numerical factorization.
  void initSolver();

  void loop();

  ///@brief set everything to starting configuration.
  void resetSim();

  void saveMeshState();

  float compressRatio;

  int compressSteps;

  ///@brief number of time steps taken.
  int     frameCnt;

  StepperStatic * stepper;

  std::string   simDirName;

  bool   savePos;
};

int loadCompressScene(WorldCompress * w, const ConfigFile & conf);

#endif