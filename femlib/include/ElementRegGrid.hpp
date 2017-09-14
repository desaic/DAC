/*
 * ElementRegGrid.hpp
 *
 *  Created on: Aug 29, 2013
 *      Author: desaic
 */

#ifndef ELEMENTREGGRID_HPP_
#define ELEMENTREGGRID_HPP_

#include "ElementMesh.hpp"
#include <Eigen/Dense>

class ElementRegGrid:public ElementMesh
{
public:
  ElementRegGrid(int _nx = 0 , int _ny = 0, int _nz = 0);
  void resize(int _nx, int _ny, int _nz);
  void allocate();
  void deleteEle();

  ///@brief remove unused vertices
  void rmEmptyVert();
  
  int GetEleInd(const Eigen::Vector3d & p);

  ///@brief 
  int GetEleInd(int ii , int jj, int kk)const;
  int GetEleInd_clamp(const Eigen::Vector3d &p);
  int GetVertInd(int ii , int jj, int kk);
  int toArrInd(int ii , int jj, int kk)const{return ii * ny * nz + jj * nz + kk;}

  virtual ~ElementRegGrid();
  std::vector<int> eleIdx;
  //map from  grid index to vertex index in array X and x.
  std::vector<int> vertIdx;
  int nx,ny,nz;
  float dx;

  Eigen::Vector3d origin;
};
#endif /* ELEMENTREGGRID_HPP_ */
