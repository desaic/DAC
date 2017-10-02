#include "CCDTrigMesh.hpp"
#include "DeformModel.h"
#include "timing.h"

void CCDTrigMesh::init(const double * verts, int nv, int * trigs, int nt)
{

  mdl = new DeformModel(verts, nv, trigs, nt, modelScale);
  TIMING_BEGIN
  mdl->Decompose();
  TIMING_END("Decompose")
  bool ccd = true;
  TIMING_BEGIN
  mdl->BuildBVH(ccd);
  TIMING_END("Build BVH")
}

void CCDTrigMesh::update(const double * verts, int nv)
{
  bool ccd = true;
  //mdl->Deform();
  mdl->copyVert(verts, nv, modelScale);
  mdl->UpdateTriNorm();
  mdl->RebuildBVH(ccd);
  mdl->ResetCounter();
  mdl->SelfCollide(ccd);
  //mdl->UpdateCollide();
}
