#ifndef ELEMENT_MESH_UTIL_HPP
#define ELEMENT_MESH_UTIL_HPP

#include "ElementMesh.hpp"

extern int HexFaces[6][4];

///@brief save abaqus format.
void saveAbq(ElementMesh * em, std::string filename, float scale = 1);

//find the vertex that would collide during loading.
//Hack to avoid collision detection
void specialVertexBound(ElementMesh * m, float lb);

void hexToTrigMesh(ElementMesh * em, TrigMesh * tm);

void hexToTrigStress(ElementMesh * em, TrigMesh * tm,
  std::vector<std::vector<double> > & P);

void updateSurfVert(const std::vector<Eigen::Vector3d> & x, TrigMesh * tm, const std::vector<int> & vidx);

std::vector<std::vector<Eigen::Vector3d> >
loadModes(const char * filename, int nvert, int dim);

void savePartObj(ElementMesh * em, int part, std::string filename);

void saveSimState(const ElementMesh * em, std::ostream & out);

void loadSimState(ElementMesh * em, std::istream & in, bool loadDisplacement = false);

/// \brief compute piola kirchoff stress at each quadrature point.
void quadratureStress(ElementMesh * em, std::vector<std::vector<double> > & P);

/// add equality constraint to stiffness matrix K.
/// @param dofs all subsequent dofs are constrained to be equal to the first dof.
/// changes matrix size by constructing kkt system.
/// remember to change x and b sizes.
void equalityConstraint(std::vector<int> dofs, ElementMesh * em,
  Eigen::SparseMatrix<double> & K);

ElementMesh* subdivideHex(ElementMesh * em, int nSub = 2);

///@brief moves rest pose X to origin
///and sets x=X.
void placeAtOrigin(ElementMesh * em);

void scaleMesh(ElementMesh * em, float scale);
///\brief scales the rest shape "X" of the mesh
void scaleMesh(ElementMesh * em, const std::vector<float> & scale);
///\brief scale the deformed mesh "x".
void scaleMeshx(ElementMesh * em, const std::vector<float> & scale);
void translate(ElementMesh * em, const Eigen::Vector3d & vec);

void savePly(TrigMesh * tm, std::string filename);

void saveVelObj(ElementMesh * em, std::string filename);
void saveStiffness(ElementMesh * em);

///@param O new corner (0,0,0) of bounding box.
///@param Frame. Each column is an axis in world space.
void applyFrame(ElementMesh * em, const Eigen::Vector3d & O,
  const Eigen::Matrix3d & Frame);

bool hitWall(float wallDist, const std::vector<Eigen::Vector3d> & x, int dim, int sign);

void saveModalMeshes(ElementMesh * em, std::string modefile);
#endif