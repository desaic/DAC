#include "EigenUtil.hpp"

Eigen::Matrix3Xd assemble3n(const std::vector<Eigen::Vector3d> & a)
{
  Eigen::Matrix3Xd M(3, a.size());
  for(unsigned int ii = 0; ii<a.size(); ii++){
    M.col(ii) = a[ii];
  }
  return M;
}

Eigen::VectorXd
assemblex(const std::vector<Eigen::Vector3d> & a, int dim)
{
  Eigen::VectorXd x(dim *a.size());
  for(int ii = 0; ii<a.size(); ii++){
    for(int jj = 0; jj<dim; jj++){
      x[ii*dim + jj] = a[ii][jj];
    }
  }
  return x;
}

void eigen2vector(const Eigen::VectorXd & ev, std::vector<double> & v)
{
  assert(ev.rows() == (int)v.size());
  for(unsigned int ii = 0; ii<v.size(); ii++){
    v[ii] = ev[ii];
  }
}

void vector2eigen(const std::vector<double> & v, Eigen::VectorXd & ev)
{
  assert(ev.rows() == (int)v.size());
  for (unsigned int ii = 0; ii<v.size(); ii++){
    ev[ii] = v[ii];
  }
}


void write_matlab(std::ostream &output, const char *variable_name,
                  const Eigen::SparseMatrix<double> & M)
{

  output<<"I=[";
   for(int i=0; i<M.cols(); ++i){
     for (Eigen::SparseMatrix<double>::InnerIterator it(M, i); it; ++it){
       output<<" "<<i+1;
     }
     if(i>0 && i % 100==0){
       output<<" ...\n";
     }
   }

   output<<"];\n  J=[";
   for(int i=0; i<M.cols(); ++i){
     for (Eigen::SparseMatrix<double>::InnerIterator it(M, i); it; ++it){
       output<<" "<<it.row()+1;
     }
     if(i>0 && i % 100==0){
       output<<" ...\n";
     }
   }

   output<<"];\n  V=[";
   for(int i=0; i<M.cols(); ++i){
     for (Eigen::SparseMatrix<double>::InnerIterator it(M, i); it; ++it){
       output<<" "<<it.value();
     }
     if(i>0 && i % 100==0){
       output<<" ...\n";
     }
   }
   output<<"];\n";
   output<<variable_name<<"=sparse(I,J,V,";
   output<<M.cols()<<", "<<M.rows()<<");"<<std::endl;
}

void write_matlab_lists(std::ostream &output,
                        const Eigen::SparseMatrix<double> & M)
{
  int nnz = M.nonZeros();
  output<<M.rows()<<" "<<nnz<<"\n";
  for(int i=0; i<M.cols(); ++i){
    for (Eigen::SparseMatrix<double>::InnerIterator it(M, i); it; ++it){
      output<< i+1 << " " << it.row()+1 << " " << it.value() <<"\n";
    }
  }
}

void write_vega_lists(std::ostream &output,
                        const Eigen::SparseMatrix<double> & M)
{
  output<<M.rows()<<"\n"<<M.cols()<<"\n";
  output.precision(16);
  for(int i=0; i<M.cols(); ++i){
    for (Eigen::SparseMatrix<double>::InnerIterator it(M, i); it; ++it){
      output<< i << " " << it.row() << " " << it.value() <<"\n";
    }
  }
}
