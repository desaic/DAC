#ifndef CCD_TRIG_MESH
#define CCD_TRIG_MESH

class DeformModel;

class CCDTrigMesh
{
public:
  CCDTrigMesh() :modelScale(1000.0f),mdl(0){}
  void init(const double * verts, int nv, int * trigs, int nt);
  void update(const double * verts, int nv);
  float modelScale;
  DeformModel * mdl;
};

#endif