#ifndef MESH_H
#define MESH_H

#include <map>
#include <vector>
#include <fstream>
#include <Eigen/Dense>

#include "MathDefines.h"
class ElementMesh;
class TrigMesh{
public:
  std::vector<Vector3s>v;
  std::vector<Vector3s>n;
  std::vector<Vector2s>tex;
  std::vector<Vector3i>texId;
  //element index corresponding to a triangle.
  std::vector<int> ei;
  //face index corresponding to a quad.
  std::vector<int> fi;
  ///@brief vertex index in element mesh corresponding to v.
  std::vector<int> vidx;
  ///@brief triangles
  std::vector<Vector3i>t;

  ///\brief triangle color for each vertex in each row.
  std::vector<Eigen::Matrix3f> tcolor;

  TrigMesh();
  TrigMesh(const std::vector<Vector3s>&_v,
      const std::vector<Vector3i>&_t);
  TrigMesh(const Vector3s *_v, const Vector3i *_t);
  
  TrigMesh(const char * filename,bool normalize);
  void load_mesh(const char * filename, bool normalize=true);
  void save(const char * filename);
  void save(std::ostream &out, std::vector<Vector3s> *vert=0);
  void load(std::istream &in);
  void read_ply(std::istream & f);
  void read_obj(std::istream &f);

  void save_obj(const char * filename);
  //void save_obj(const char * filename, int mat, ElementMesh * em);
  void load_tex(const char * filename);
  
  void compute_norm();
  void rescale();
  void append(const TrigMesh & m);
  TrigMesh & operator= (const TrigMesh& m);
  virtual void update();
};


void makeCube(TrigMesh & m, const Vector3s &mn,
    const Vector3s mx);
///@brief cube [0,1]^3
extern TrigMesh UNIT_CUBE;

void adjlist(const TrigMesh & m, std::vector<std::vector<int> > & adjMat);

#endif
