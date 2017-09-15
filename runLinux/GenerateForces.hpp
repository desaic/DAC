#ifndef GENERATE_FORCES_HPP
#define GENERATE_FORCES_HPP
#include "Eigen/Dense"
#include <vector>

class ForceSample
{
public:
  ///@brief Forces are ordered according to ElementHex
  ForceSample();
  ///@brief generate force bases that can be combined later
  void generateBase();
  ///@brief generate force for each vertex in each axis
  void generateAxisForce();
  ///@brief generate force of other kinds
  void generateModeForce();
  void generateTwist(bool flip = false);

  ///@brief add positive and negative base forces
  void addBases();

  void saveForces(const char * filename);
  void loadForces(const char * filename);
  std::vector<std::vector<Eigen::Vector3d> > baseForce;
  std::vector<std::vector<Eigen::Vector3d> > ff;
  void scaleForces();
  float forceScale;
  //number of steps in each force direction
  int nMag = 5;
};

#endif