#include "ElementMesh.hpp"
#include "Element.hpp"
#include "Eigen/Sparse"
#include "Timer.hpp"

#include <iostream>
#include <fstream>

typedef Eigen::Triplet<double> Tripletd;

Eigen::SparseMatrix<double>
ElementMesh::getStiffnessSparse()
{
  int N = dim * (int)x.size();
  Timer timer;
  if (Kpattern.size() == 0){
    std::vector<int> I, J;
    stiffnessPattern(I, J);
  }

  int nV = e[0]->nV();
  if (sblocks.size() == 0){
    sblocks.resize(x.size());
    for (unsigned int ii = 0; ii < e.size(); ii++){
      Element * ele = e[ii];
      nV = ele->nV();
      for (int jj = 0; jj < nV; jj++){
        int vj = ele->at(jj);
        for (int kk = 0; kk < nV; kk++){
          int vk = ele->at(kk);
          auto it = sblocks[vj].find(vk);
          if (it == sblocks[vj].end()){
            sblocks[vj][vk] = Eigen::MatrixXd::Zero(dim, dim);
          }
        }
      }
    }
  }
  for (unsigned int vi = 0; vi < x.size(); vi++){
    for (auto it = sblocks[vi].begin(); it != sblocks[vi].end(); it++){
      it->second = Eigen::MatrixXd::Zero(dim,dim);
    }
  }

  double localStiffTime = 0;
  double assembleTime = 0;
  for(unsigned int ii = 0;ii<e.size();ii++){
    timer.startWall();
    Element * ele = e[ii];
    Eigen::MatrixXd K  = getStiffness(ii);
    timer.endWall();
    localStiffTime += timer.getSecondsWall();
    timer.startWall();
    for (int jj = 0; jj<ele->nV(); jj++){
      int vj = ele->at(jj);
      for (int kk = 0; kk<ele->nV(); kk++){
        int vk = ele->at(kk);
        auto it = sblocks[vj].find(vk);
        it->second += K.block(dim * jj, dim * kk, dim, dim);
      }
    }
    timer.endWall();
    assembleTime += timer.getSecondsWall();
  }
  //std::cout << "local stiffness time " << localStiffTime << "\n";
  //std::cout << "assemble time " << assembleTime << "\n";
  timer.startWall();
  int cnt = 0;
  for(unsigned int vi = 0; vi<x.size(); vi++){
    for(int di = 0; di<dim; di++){
      int col = vi*dim + di;
      for (auto it = sblocks[vi].begin(); it != sblocks[vi].end(); it++){
        int vj = it->first;
        for(int dj = 0; dj<dim; dj++){
          int row = vj*dim+dj;
#ifdef DEBUG
          if(row>=N || col>=N){
            std::cout<<"Index out of bound in stiffness assembly" <<col<<" "<<row<<"\n";
          }
#endif
          Kpattern.valuePtr()[cnt] = (it->second)(di, dj);
          cnt++;
        }
      }
    }
  }

  timer.endWall();
  //std::cout<<"Stiffness conversion time "<<timer.getSecondsWall()<<"\n";
  return Kpattern;
}

Eigen::SparseMatrix<double>
ElementMesh::getStiffnessDamping(Eigen::SparseMatrix<double> & D,
  std::vector<float> dcoef)
{
  int N = dim * (int)x.size();
  Timer timer;
  if (Kpattern.size() == 0){
    std::vector<int> I, J;
    stiffnessPattern(I, J);
  }
  D = Kpattern;
  int nV = e[0]->nV();
  if (sblocks.size() == 0){
    sblocks.resize(x.size());
    for (unsigned int ii = 0; ii < e.size(); ii++){
      Element * ele = e[ii];
      nV = ele->nV();
      for (int jj = 0; jj < nV; jj++){
        int vj = ele->at(jj);
        for (int kk = 0; kk < nV; kk++){
          int vk = ele->at(kk);
          auto it = sblocks[vj].find(vk);
          if (it == sblocks[vj].end()){
            sblocks[vj][vk] = Eigen::MatrixXd::Zero(dim, dim);
          }
        }
      }
    }
  }
  for (unsigned int vi = 0; vi < x.size(); vi++){
    for (auto it = sblocks[vi].begin(); it != sblocks[vi].end(); it++){
      it->second = Eigen::MatrixXd::Zero(dim, dim);
    }
  }
  std::vector< std::map< int, Eigen::MatrixXd > > Dblocks = sblocks;

  double localStiffTime = 0;
  double assembleTime = 0;

  for (unsigned int ii = 0; ii<e.size(); ii++){
    timer.startWall();
    int mat = me[ii];
    float damp = dcoef[mat];
    Element * ele = e[ii];
    nV = ele->nV();
    Eigen::MatrixXd K = getStiffness(ii);
    timer.endWall();
    localStiffTime += timer.getSecondsWall();
    timer.startWall();
    for (int jj = 0; jj<nV; jj++){
      int vj = ele->at(jj);
      for (int kk = 0; kk<nV; kk++){
        int vk = ele->at(kk);
        auto it = sblocks[vj].find(vk);
        it->second += K.block(dim * jj, dim * kk, dim, dim);
        it = Dblocks[vj].find(vk);
        it->second += damp * K.block(dim * jj, dim * kk, dim, dim);
      }
    }
    timer.endWall();
    assembleTime += timer.getSecondsWall();
  }
  //std::cout << "local stiffness time " << localStiffTime << "\n";
  //std::cout << "assemble time " << assembleTime << "\n";
  timer.startWall();
  int cnt = 0;
  for (unsigned int vi = 0; vi<x.size(); vi++){
    for (int di = 0; di<dim; di++){
      int cnt1 = cnt;
      for (auto it = sblocks[vi].begin(); it != sblocks[vi].end(); it++){
        int vj = it->first;
        for (int dj = 0; dj<dim; dj++){
          int row = vj*dim + dj;
          Kpattern.valuePtr()[cnt1] = (it->second)(di, dj);
          cnt1++;
        }
      }

      for (auto it = Dblocks[vi].begin(); it != Dblocks[vi].end(); it++){
        int vj = it->first;
        for (int dj = 0; dj<dim; dj++){
          int row = vj*dim + dj;
          D.valuePtr()[cnt] = (it->second)(di, dj);
          cnt++;
        }
      }

    }
  }

  timer.endWall();
  //std::cout<<"Stiffness conversion time "<<timer.getSecondsWall()<<"\n";
  return Kpattern;

}

void ElementMesh::stiffnessPattern(std::vector<int> & I, std::vector<int> & J)
{
  int N = dim * (int)x.size();
  Kpattern.resize(N, N);
  std::cout << "N " << N << "\n";
  Kpattern.reserve(Eigen::VectorXi::Constant(N, 81));
  for(unsigned int ii = 0;ii<e.size();ii++){
    Element * ele = e[ii];
    for(int jj = 0; jj<ele->nV(); jj++){
      int vj = ele->at(jj);
      for(int kk = 0; kk<ele->nV(); kk++){
        int vk = ele->at(kk);
        for(int dim1= 0 ;dim1<dim; dim1++){
          for(int dim2= 0 ;dim2<dim; dim2++){
            Kpattern.coeffRef(dim*vj + dim1, dim*vk + dim2) += 1;
          }
        }
      }
    }
  }
  Kpattern.makeCompressed();
  I.push_back(0);
  for (int ii = 0; ii<Kpattern.cols(); ii++){
    for (Eigen::SparseMatrix<double>::InnerIterator it(Kpattern, ii); it; ++it){
      J.push_back(it.row());
      //std::cout << it.row() << "\n";
    }
   I.push_back((int)J.size());
  }
}

void fixTranslation(Eigen::SparseMatrix<double> & K, ElementMesh * mesh)
{
  int nrow = K.rows();
  int ncol = K.cols();
  int nconstraints = 3;

  K.conservativeResize(nrow + nconstraints, ncol + nconstraints);
  Eigen::SparseMatrix<double> KConstraint(nrow + nconstraints, nrow + nconstraints);

  Eigen::Vector3d center(0, 0, 0);
  for (int ii = 0; ii<mesh->x.size(); ii++){
    center += mesh->x[ii];
  }
  center /= (double)mesh->x.size();
  double cScale = 1.0;
  std::vector<Tripletd> triplets;
  for (int ii = 0; ii<mesh->x.size(); ii++)
  {
    int col = 3 * ii;
    //relative position to center
    Eigen::Vector3d rel = mesh->x[ii] - center;

    //fix translation
    for (int kk = 0; kk<3; kk++){
      triplets.push_back(Tripletd(nrow + kk, col + kk, -cScale));
      triplets.push_back(Tripletd(col + kk, nrow + kk, cScale));
    }
  }
  for (int ii = 0; ii<nconstraints; ii++){
    triplets.push_back(Tripletd(nrow + ii, nrow + ii, 0.f));
  }
  KConstraint.setFromTriplets(triplets.begin(), triplets.end());
  K += KConstraint;
}

void fixRotation(Eigen::SparseMatrix<double> & K, ElementMesh * mesh)
{
  int nrow = K.rows();
  int ncol = K.cols();
  int nconstraints = 3;

  K.conservativeResize(nrow + nconstraints, ncol + nconstraints);
  Eigen::SparseMatrix<double> KConstraint(nrow + nconstraints, nrow + nconstraints);

  Eigen::Vector3d center(0, 0, 0);
  for (int ii = 0; ii<mesh->x.size(); ii++){
    center += mesh->x[ii];
  }
  center /= (double)mesh->x.size();
  double cScale = 1.0;
  std::vector<Tripletd> triplets;
  for (int ii = 0; ii<mesh->x.size(); ii++){
    int col = 3 * ii;
    //relative position to center
    Eigen::Vector3d rel = mesh->x[ii] - center;
    Eigen::Matrix3d c;
    c << 0, -rel[2], rel[1],
      rel[2], 0, -rel[0],
      -rel[1], rel[0], 0;
    c = -cScale * c;

    // fix rotation
    for (int kk = 0; kk<3; kk++){
      for (int ll = 0; ll<3; ll++){
        triplets.push_back(Tripletd(nrow + kk, col + ll, c(kk, ll)));
        triplets.push_back(Tripletd(col + ll, nrow + kk, -c(kk, ll)));
      }
    }
  }
  for (int ii = 0; ii<nconstraints; ii++){
    triplets.push_back(Tripletd(nrow + ii, nrow + ii, 0.f));
  }
  KConstraint.setFromTriplets(triplets.begin(), triplets.end());
  K += KConstraint;
}
