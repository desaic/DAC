// DMVMap.h
//
// Breannan Smith
// Last updated: 09/22/2015

#ifndef DMV_MAP_H
#define DMV_MAP_H
#include "MathDefines.h"
class RigidBody3DSystem;
class DMVMap
{

public:

  virtual ~DMVMap();

  virtual void flow( const VectorXs& q0, const VectorXs& v0, RigidBody3DSystem& fsys, const unsigned iteration, const scalar& dt, VectorXs& q1, VectorXs& v1 );

  virtual std::string name() const;

};

#endif
