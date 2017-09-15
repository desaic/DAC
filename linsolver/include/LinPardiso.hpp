#ifndef LINPARDISO_HPP
#define LINPARDISO_HPP

#include "LinSolver.hpp"
#include "pardiso_sym.hpp"

template <typename T>
class LinPardiso:public LinSolver<T>
{

public:

  std::vector<int> ia;
  std::vector<int> ja;
  std::vector<double> val;
  bool triangle;

  PardisoState state;

  ~LinPardiso(){
    if(ia.size()>0){
      free();
    }
  }

  virtual void init(){
    state.mtype=-2;
    pardisoInit(&state);
    state.maxfct = 1;		/* Maximum number of numerical factorizations.  */
    state.mnum   = 1;         /* Which factorization to use. */
    state.msglvl = 0;//1;         /* Print statistical information  */
    //only store triangle for symmetric system
    triangle = true;
  }

  ///@brief symbolic factorization.
  void init(const Eigen::SparseMatrix<T> & A){
    //have been intialized before.
    if(ia.size()>0){
      free();
      pardisoInit(&state);
    }
    int ncol = A.cols();
    ia.resize(ncol+1,0);
    ja.resize(A.nonZeros());
    int cnt = 0;
    for(int ii = 0; ii<ncol; ii++){
      for (Eigen::SparseMatrix<double>::InnerIterator it(A, ii); it; ++it){
        //only add lower triangle
        int row = it.row();
        if(triangle && ii>row){
          continue;
        }
        ja[cnt] = row;
        cnt++;
      }
      ia[ii+1] = cnt;
    }

    /* -------------------------------------------------------------------- */
    /* ..  Convert matrix from 0-based C-notation to Fortran 1-based        */
    /*     notation.                                                        */
    /* -------------------------------------------------------------------- */
    for (unsigned int i = 0; i < ia.size(); i++) {
      ia[i] += 1;
    }
    for (unsigned int i = 0; i < ja.size(); i++) {
      ja[i] += 1;
    }
    pardisoSymbolicFactorize(&ia[0], &ja[0], ncol, &state);
  }

  void free(){
    int ncol = (int)ia.size()-1;
    pardisoFree(&ia[0], &ja[0], ncol, &state);
    ia.clear();
    ja.clear();
    val.clear();
  }

  ///@return status. >= 0 For success. <0 for error codes.
  virtual int solve(Eigen::SparseMatrix<T> & A,
                    const std::vector<T> & b, std::vector<T> & x,
                    T & residual, int iters){
    int ncol = (int)ia.size()-1;
    val.resize(A.nonZeros());
    //cast T to double for pardiso
    std::vector<double> rhs = b;
    int valIdx = 0;
    for(int ii = 0; ii<ncol; ii++){
      for (Eigen::SparseMatrix<double>::InnerIterator it(A, ii); it; ++it){
        //only add lower triangle
        int row = it.row();
        if(triangle && ii>row){
          continue;
        }
        val[valIdx]=it.value();
        valIdx++;
      }
    }

    int status = pardisoSolve( &(ia[0]), &(ja[0]), &(val[0]), ncol,
        &(x[0]), &(rhs[0]), &state);

    Eigen::Map< Eigen::Matrix<T,Eigen::Dynamic,1> >ev(x.data(), x.size());
    Eigen::Matrix<T,Eigen::Dynamic,1> Ax = A*ev;

    residual =0;
    double r0 = 1e-20;
    for(unsigned int ii =0 ; ii<b.size(); ii++){
      double diff = Ax[ii] - b[ii];
//      std::cout<<Ax[ii]<<" "<<b[ii]<<"\n";
      residual += diff*diff;
      r0 += b[ii] * b[ii];
    }
    residual = std::sqrt(residual/r0);
    iters= 0;
    return status;
  }

  ///@brief back substitute using precomputed factorization.
  virtual int solve(T * b, std::vector<T> &x, T &residual, int iters){
    int ncol = (int)x.size();

    //make sure rhs has the same size the matrix that has been factorized.
    assert(ncol == (int)ia.size()-1);

    int status = pardisoBackSubstitute( &ia[0],&ja[0],&val[0], ncol,
        &(x[0]), b, &state);
    return status;
  }

};

#endif
