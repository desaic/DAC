#include "Element.hpp"
#include "ElementMesh.hpp"
#include "Material.hpp"


Eigen::MatrixXd
ElementMesh::getStiffness(int eIdx)
{
  Eigen::MatrixXd K = m[me[eIdx]]->getStiffness(eIdx,this);
  return K;
}

Eigen::MatrixXd
ElementMesh::getStiffness()
{
  int matSize = dim * (int)x.size();
  Eigen::MatrixXd Kglobal=Eigen::MatrixXd::Zero(matSize,matSize);
  for(unsigned int ii = 0;ii<e.size();ii++){
    Eigen::MatrixXd K = getStiffness(ii);
    for(int jj = 0; jj<e[ii]->nV(); jj++){
      int vj = e[ii]->at(jj);
      for(int kk = 0; kk<e[ii]->nV(); kk++){
        int vk = e[ii]->at(kk);    
        Kglobal.block(dim*vj, dim*vk, dim, dim) +=
            K.block(dim*jj, dim*kk, dim, dim);
      }
    }
  }
  return Kglobal;
}
