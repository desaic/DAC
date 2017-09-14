#ifndef EIGENUTIL_HPP
#define EIGENUTIL_HPP

#include <Eigen/Sparse>
#include <Eigen/SVD>
#include <vector>
#include <iostream>

Eigen::Matrix3Xd assemble3n(const std::vector<Eigen::Vector3d> & a);

Eigen::VectorXd assemblex(const std::vector<Eigen::Vector3d> & a, int dim = 2);

void eigen2vector(const Eigen::VectorXd & ev, std::vector<double> & v);
void vector2eigen(const std::vector<double> & v, Eigen::VectorXd & ev);
///@brief write sparse matrix in matlab format such that one can copy paste
///into matlab workspace
void write_matlab(std::ostream &output, const char *variable_name,
                  const Eigen::SparseMatrix<double> & M);

///@brief write sparse matrix in matlab format for loading with matlab code.
void write_matlab_lists(std::ostream &output,
                        const Eigen::SparseMatrix<double> & M);

///@brief write sparse matrix in Vega format for loading with vega sparsematrix.
void write_vega_lists(std::ostream &output,
                        const Eigen::SparseMatrix<double> & M);


///@brief zero out off-diagonal terms of row and col ii if fixed[ii]!=0.
template <typename T>
void zeroOffDiag(Eigen::SparseMatrix<T> & K, const std::vector<int> & fixed)
{
  for(int col =0; col<K.cols(); col++){
    for (Eigen::SparseMatrix<double>::InnerIterator it(K, col); it; ++it){
      int row = it.row();
      if( (fixed[col] != 0 || fixed[row]!=0) && row != col){
        it.valueRef() = 0;
      }
    }
  }
}

///@param x 3xn matrix of deformed verts
///@param p 3xn matrix of ref verts
template <typename T>
Eigen::Matrix<T,3,3> shapeMatch(const Eigen::Matrix<T, 3, Eigen::Dynamic> & x,
                           const Eigen::Matrix<T, 3, Eigen::Dynamic> & p)
{
  Eigen::Matrix<T, 3, 1> xcm, pcm;

  Eigen::Matrix<T, 3, Eigen::Dynamic> Mx, Mp;

  assert(x.cols() == p.cols());
  int cols = x.cols();
  for(int row = 0; row<3; row++){
    xcm[row] = x.row(row).sum()/x.cols();
    pcm[row] = p.row(row).sum()/p.cols();
  }
  Mx = x-xcm.replicate(1, cols);
  Mp = p-pcm.replicate(1, cols);

  Eigen::Matrix<T, 3, 3> A, A1, A2;
  A1 = Mp*Mx.transpose();
  A2 = Mx*Mx.transpose();
  A = A1 * A2.inverse();
  Eigen::JacobiSVD<Eigen::Matrix<T,3,3> > svdA(A, Eigen::ComputeFullU | Eigen::ComputeFullV);

  Eigen::Matrix<T,3,3> Vt = svdA.matrixV().transpose();
  Eigen::Matrix<T,3,3> R = svdA.matrixU()*Vt;
  return R;
}

template <typename T>
Eigen::Matrix<T,3,3> shapeMatch2d(const Eigen::Matrix<T, 3, Eigen::Dynamic> & x,
                           const Eigen::Matrix<T, 3, Eigen::Dynamic> & p)
{
  Eigen::Matrix<T, 3, 1> xcm, pcm;

  Eigen::Matrix<T, 3, Eigen::Dynamic> Mx, Mp;

  assert(x.cols() == p.cols());
  int cols = x.cols();
  for(int row = 0; row<3; row++){
    xcm[row] = x.row(row).sum()/x.cols();
    pcm[row] = p.row(row).sum()/p.cols();
  }
  Mx = x-xcm.replicate(1, cols);
  Mp = p-pcm.replicate(1, cols);

  Eigen::Matrix<T, 3, 3> A, A1, A2;

  A1 = Mp*Mx.transpose();
  A2 = Mx*Mx.transpose();
  A1(2,2) = 1;
  A2(2,2) = 1;
  A = A1 * A2.inverse();
  Eigen::JacobiSVD<Eigen::Matrix<T,3,3> > svdA(A, Eigen::ComputeFullU | Eigen::ComputeFullV);

  Eigen::Matrix<T,3,3> Vt = svdA.matrixV().transpose();
  Eigen::Matrix<T,3,3> R = svdA.matrixU()*Vt;
  return R;
}

template <typename T>
Eigen::SparseMatrix<T>
rmConstrained(const Eigen::SparseMatrix<T> & K,
              const std::vector<int> & fixed)
{
  int cnt = 0;
  std::vector<int> cidx(K.rows(),0);
  for(unsigned int ii = 0; ii<fixed.size(); ii++){
    if(fixed[ii]){
      continue;
    }
    cidx[ii] = cnt;
    cnt++;
  }

  int N = cnt;
//  std::cout<<N<<"\n";
  Eigen::SparseMatrix<T> Kc(N, N);

  int max_nnz = 0;
  for(unsigned int col =0; col<K.cols(); col++){
    int n = 0;
    for (Eigen::SparseMatrix<double>::InnerIterator it(K, col); it; ++it){
      n++;
    }
    if(n>max_nnz){
      max_nnz = n;
    }
  }
  Kc.reserve(Eigen::VectorXi::Constant(N, max_nnz));

  for(unsigned int col =0; col<K.cols(); col++){
    for (Eigen::SparseMatrix<double>::InnerIterator it(K, col); it; ++it){
      int row = it.row();
      if(fixed[col] || fixed[row]){
        continue;
      }
      Kc.insert(cidx[row], cidx[col]) = it.value();
    }
  }

  return Kc;
}


#endif // EIGENUTIL_HPP
