// MomentTools.h
//
// Breannan Smith
// Last updated: 09/14/2015

#ifndef MOMENT_TOOLS_H
#define MOMENT_TOOLS_H

#include "MathDefines.h"

class ElementMesh;

namespace MomentTools
{

  // Implementation of Polyhedral Mass Properties (Revisited) by David Eberly
  // *** Caller must rescale mass and I by density ***
  void computeMoments( const Matrix3Xsc& vertices, const Matrix3Xuc& indices, scalar& mass, Vector3s& I, Vector3s& center, Matrix3s& R );
  void computeMoments( ElementMesh * em, Vector3s& I, Vector3s& center, Matrix3s& R );
}
void diagonalizeInertiaTensor( const Matrix3s& I, Matrix3s& R0, Vector3s& I0 );
#endif
