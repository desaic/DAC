#ifndef LINSOLVER_HPP
#define LINSOLVER_HPP

#include <Eigen/Sparse>
#include <vector>
///@brief wrapper class for calling external linear solvers on
///eigen sparse matrix (float or double precision).

template <typename T>
class LinSolver{

public:

  virtual void init(){}

  ///@brief initialization that depends on the matrix
  virtual void init(const Eigen::SparseMatrix<T> & A){}

  virtual ~LinSolver(){}

  ///@return status. >= 0 For success. <0 for error codes.
  virtual int solve(Eigen::SparseMatrix<T> & A,
                    const std::vector<T> & b, std::vector<T> & x,
                    T & residual, int iters){return -1;}

  virtual int compute(Eigen::SparseMatrix<T> & A){ 
    //default unimplemented
    return -1; }

  ///@brief solve Ax=b where A has already been computed e.g. numerical factorized by the solver.
  virtual int solve(T * b, std::vector<T> &x, T &residual, int iters){
    //default unimplemented
    return -1;
  }

};

#endif
