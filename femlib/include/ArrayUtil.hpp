#ifndef ARRAYUTIL_HPP
#define ARRAYUTIL_HPP

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <Eigen/Dense>
///@brief for Vector3f/d
template <typename T>
int
save_arr(const std::vector<T> & a, std::ostream & out)
{
  out<<a.size()<<"\n";
  for(unsigned int ii = 0; ii<a.size(); ii++){
    for(int jj = 0; jj<a[ii].size(); jj++){
      if(jj>0){
        out<<" ";
      }
      out<<a[ii][jj];
    }
    out<<"\n";
  }
  return 0;
}

int
load_arr(std::vector<Eigen::Vector3d> & a, std::istream & in);

///@param a bilinear weights in [0 1].
void bilinearWeights(float * a, float * w);

template <typename T>
void printVector(const std::vector<T> & v, std::ostream & out)
{
  for(unsigned int ii =0 ; ii<v.size(); ii++){
    out<<v[ii]<<"\n";
  }
}

template <typename T>
void BBox(const std::vector<T>& v, T & mn, T & mx)
{
  mn = v[0];
  mx = v[0];
  int DIM = (int)v[0].size();
  for(unsigned int ii = 1 ;ii<v.size();ii++){
    for(int dim = 0; dim < DIM; dim++){
      if(v[ii][dim]<mn[dim]){
        mn[dim] = v[ii][dim];
      }
      if(v[ii][dim]>mx[dim]){
        mx[dim] = v[ii][dim];
      }
    }
  }
}

template<typename T>
void add(std::vector<T> & dst, const std::vector<T> & src)
{
  for(unsigned int ii = 0;ii<dst.size();ii++){
    dst[ii] += src[ii];
  }
}

void add_scaled(std::vector<double> & dst, float scale,
  const std::vector<double> & src);

double infNorm(const std::vector<double> & a);

template<typename T>
std::vector<T> mul(float f, const std::vector<T> & src)
{
  std::vector<T> prod(src.size());
  for(unsigned int ii = 0;ii<src.size();ii++){
    prod[ii] = f*src[ii];
  }
  return prod;
}

template<typename T>
void addmul(std::vector<T> & dst, float f, const std::vector<T> & src)
{
  for(unsigned int ii = 0;ii<src.size();ii++){
    dst[ii] += f*src[ii];
  }
}

template<typename T>
void setVal(std::vector<T> & a, const std::vector<int> & c, const T & val)
{
  for(unsigned int ii =0 ; ii<a.size(); ii++){
    if(c[ii]){
      a[ii] = val;
    }
  }
}

int gridToLinearIdx(int ix, int iy, int iz, const std::vector<int> & gridSize);

int find(int val, const std::vector<int> & a);
#endif // ARRAYUTIL_HPP
