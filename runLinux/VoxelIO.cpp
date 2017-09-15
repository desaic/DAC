#include "VoxelIO.hpp"
#include "ArrayUtil.hpp"
#include "ElementRegGrid.hpp"
#include "Element.hpp"
#include "FileUtil.hpp"
#include <vector>
#include <string>
#include <iostream>
#include <fstream>

int loadArr3dTxt(std::string filename, std::vector<double> & arr,
  int &inputSize) {
  int N = 0;
  FileUtilIn in(filename.c_str());
  if (!in.good()) {
    return -1;
  }
  in.in >> N >> N >> N;
  inputSize = N;
  int nVox = N*N*N;
  const int MAX_BUF = 64;
  arr.resize(nVox, 0);
  std::cout << "loadArr3dTxt " << N << "\n";
  for (int i = 0; i < nVox; i++) {
    char inStr[MAX_BUF];
    double val = 0;
    in.in >> inStr;
    val = atof(inStr);
    arr[i] = val;
  }
  in.close();
  return 0;
}

void printIntStructure(const double * s, const std::vector<int> & gridSize,
  std::ostream & out)
{
  if (gridSize.size() < 3) {
    return;
  }
  int N = gridSize[0] * gridSize[1] * gridSize[2];
  out << gridSize[0] << " " << gridSize[1] << " " << gridSize[2] << "\n";
  for (int i = 0; i < N; i++) {
    int val = 1;
    if (s[i] < 0.5) {
      val = 0;
    }
    out << val << " ";
    if (i % gridSize[1] == gridSize[1] - 1) {
      out << "\n";
    }
  }
  out << "\n";
}

void loadBinaryStructure(const std::string & filename,
  std::vector<int> & s,
  std::vector<int> & gridSize)
{
  int dim = 3;
  std::ifstream in(filename, std::ios::in | std::ios::binary);
  if (!in.good()) {
    std::cout << "Cannot open input " << filename << "\n";
    in.close();
    return;
  }
  gridSize.clear();
  gridSize.resize(dim, 0);
  int nCell = 1;
  for (int i = 0; i < dim; i++) {
    in.read((char*)(&gridSize[i]), sizeof(int));
    nCell *= gridSize[i];
  }
  s.resize(nCell, 0);
  for (std::vector<bool>::size_type i = 0; i < nCell;) {
    unsigned char aggr;
    in.read((char*)&aggr, sizeof(unsigned char));
    for (unsigned char mask = 1; mask > 0 && i < nCell; ++i, mask <<= 1)
      s[i] = (aggr & mask)>0;
  }
  in.close();
}

void loadBinaryStructure(const std::string & filename,
  std::vector<double> & s,
  std::vector<int> & gridSize)
{
  std::vector<int> a;
  loadBinaryStructure(filename, a, gridSize);
  s.resize(a.size(), 0.0);
  for (size_t i = 0; i < a.size(); i++) {
    s[i] = a[i];
  }
}


std::vector<double> mirrorOrthoStructure(const std::vector<double> &s, std::vector<int> & gridSize)
{
  std::vector<int> newSize = gridSize;
  int nEle = 1;
  for (size_t i = 0; i < newSize.size(); i++) {
    newSize[i] *= 2;
    nEle *= newSize[i];
  }
  std::vector<double> t(nEle);
  for (int i = 0; i < newSize[0]; i++) {
    for (int j = 0; j < newSize[1]; j++) {
      for (int k = 0; k < newSize[2]; k++) {
        int newIdx = gridToLinearIdx(i, j, k, newSize);
        int i0 = i;
        int j0 = j;
        int k0 = k;
        if (i0 >= gridSize[0]) {
          i0 = 2 * gridSize[0] - i0 - 1;
        }
        if (j0 >= gridSize[1]) {
          j0 = 2 * gridSize[1] - j0 - 1;
        }
        if (k0 >= gridSize[2]) {
          k0 = 2 * gridSize[2] - k0 - 1;
        }
        int oldIdx = gridToLinearIdx(i0, j0, k0, gridSize);
        t[newIdx] = s[oldIdx];
      }
    }
  }
  gridSize = newSize;
  return t;
}

void getCubicTet(std::vector<double>& s, std::vector<int>& gridSize)
{
  for (int i = 0; i < gridSize[0]; i++) {
    for (int j = 0; j < gridSize[1]; j++) {
      for (int k = 0; k < gridSize[2]; k++) {
        int lidx = gridToLinearIdx(i, j, k, gridSize);
        //if (i >= gridSize[0] / 2 || j >= gridSize[1] / 2 || k >= gridSize[2] / 2) {
        //  s[lidx] = 0;
        //}
        if (j>i + 1 || k>j + 1) {
          s[lidx] = 0;
        }
      }
    }
  }
}

void loadBinDouble(const std::string & filename,
  std::vector<double>&s,
  std::vector<int> & gridSize)
{
  int dim = 3;
  std::ifstream in(filename, std::ios::in | std::ios::binary);
  if (!in.good()) {
    std::cout << "Cannot open input " << filename << "\n";
    in.close();
    return;
  }
  gridSize.clear();
  gridSize.resize(dim, 0);
  int nCell = 1;
  for (int i = 0; i < dim; i++) {
    in.read((char*)(&gridSize[i]), sizeof(int));
    nCell *= gridSize[i];
  }
  s.resize(nCell, 0);
  std::vector<double> s0(nCell);
  for (int i = 0; i < nCell; i++) {
    float val;
    in.read((char*)&val, sizeof(float));
    s0[i] = val;
  }
  in.close();

  s = mirrorOrthoStructure(s0, gridSize);

}

void saveBinaryStructure(const std::string & filename,
  const std::vector<int> & s,
  const std::vector<int> & gridSize)
{
  std::ofstream out(filename, std::ios::out | std::ios::binary);
  if (!out.good()) {
    std::cout << "Cannot open output " << filename << "\n";
    out.close();
    return;
  }
  for (size_t i = 0; i < gridSize.size(); i++) {
    out.write((const char *)&gridSize[i], sizeof(int));
  }
  for (size_t i = 0; i < s.size();) {
    unsigned char aggr = 0;
    for (unsigned char mask = 1; mask > 0 && i < s.size(); ++i, mask <<= 1)
      if (s[i])
        aggr |= mask;
    out.write((const char*)&aggr, sizeof(unsigned char));
  }
  out.close();
}

void saveBinaryStructure(const std::string & filename,
  const std::vector<double> & s,
  const std::vector<int> & gridSize)
{
  std::vector<int> a(s.size());
  float thresh = 0.5f;

  for (size_t i = 0; i < a.size(); i++) {
    a[i] = s[i] > thresh;
  }
  saveBinaryStructure(filename, a, gridSize);
}


void voxel2em(std::string inname, std::string outname)
{
  std::vector<double> s;

  std::vector<int> inputSize;

  loadBinaryStructure(inname, s, inputSize);
  std::vector<int> gridSize = inputSize;
  //for (int i = 0; i < (int)inputSize.size(); i++){
  //  gridSize[i] = 2 * inputSize[i];
  //}
  //s = upsampleVector(s, inputSize[0], inputSize[1], inputSize[2],
  //  gridSize[0], gridSize[1], gridSize[2]);
  //s = smooth3D(s, gridSize, 3);
//  s = smooth3D(s, gridSize, 3);
  
  //s = mirrorOrthoStructure(s, gridSize);
  float matThresh = 0.5f;
  ElementRegGrid * em = new ElementRegGrid();
  em->nx = gridSize[0];
  em->ny = gridSize[1];
  em->nz = gridSize[2];
  em->allocate();
  std::vector<Element*> newEle;
  int idx = 0;
  for (int i = 0; i < em->nx; i++){
    for (int j = 0; j < em->ny; j++){
      for (int k = 0; k < em->nz; k++){
        int lidx = gridToLinearIdx(i, j, k, gridSize);
        double mat = s[lidx];
        if (mat < matThresh){
          delete em->e[idx];
        }
        else{
          newEle.push_back(em->e[idx]);
        }
        idx++;
      }
    }
  }
  em->e = newEle;
  em->rmEmptyVert();
  FileUtilOut out(outname);
  em->saveMesh(out.out);
  out.close();
}

bool inbound(int i, int j, int k, const std::vector<int> & s)
{
  return (i >= 0 && i < s[0] && j >= 0 && j < s[1] &&
    k >= 0 && k < s[2]);
}

void linearToGridIdx(int l, const std::vector<int> & gridSize,
  std::vector<int> & gridIdx)
{
  gridIdx.resize(gridSize.size(), 0);
  for (int i = (int)gridSize.size() - 1; i > 0; i--){
    gridIdx[i] = l % (gridSize[i]);
    l = l / (gridSize[i]);
  }
  gridIdx[0] = l;
}

std::vector<double> smooth3D(const std::vector<double> & x, const std::vector<int> & gridSize,
  double R)
{
  std::vector<double> newx(x.size());
  int dim = 3;
  if ((int)gridSize.size() < dim || (x.size() != gridSize[0] * gridSize[1] * gridSize[2])){
    std::cout << "Smooth 3D array wrong size.\n";
  }
  int nbrSize = (int)(R);
  if (nbrSize < 1){
    nbrSize = 1;
  }
  std::vector<int> gridIdx(3);
  for (int ix = 0; ix < x.size(); ix++){
    linearToGridIdx(ix, gridSize, gridIdx);
    double val = 0;
    double weightSum = 0;
    //std::cout << gridIdx[0] << " " << gridIdx[1] << " " << gridIdx[2] << "\n";
    for (int i = gridIdx[0] - nbrSize; i <= gridIdx[0] + nbrSize; i++){
      for (int j = gridIdx[1] - nbrSize; j <= gridIdx[1] + nbrSize; j++){
        for (int k = gridIdx[2] - nbrSize; k <= gridIdx[2] + nbrSize; k++){
          if (!inbound(i, j, k, gridSize)){
            continue;
          }
          int nbrLinearIdx = gridToLinearIdx(i, j, k, gridSize);
          int nbrGridIdx[3] = { i, j, k };
          double dist = 0;
          for (int d = 0; d < dim; d++){
            double diff = nbrGridIdx[d] - gridIdx[d];
            dist += diff * diff;
          }
          dist = std::sqrt(dist);
          double w = std::pow(R - dist, 2);
          if (dist > R){
            continue;
          }
          val += w * x[nbrLinearIdx];
          weightSum += w;
        }
      }
    }
    if (x[ix] < 0.9){
      newx[ix] = val / weightSum;
    }
    else{
      newx[ix] = x[ix];
    }
  }
  return newx;
}

///@param nx size of input 3D array
std::vector<double> upsampleVector(const std::vector<double> & a, int nx, int ny, int nz,
  int nx1, int ny1, int nz1)
{
  std::vector<double> x(nx1 * ny1 * nz1);
  for (int i = 0; i < nx1; i++){
    for (int j = 0; j < ny1; j++){
      for (int k = 0; k < nz1; k++){
        int i0 = (int)((i + 0.5) / nx1*nx);
        int j0 = (int)((j + 0.5) / ny1*ny);
        int k0 = (int)((k + 0.5) / nz1*nz);
        x[i * ny1 * nz1 + j * nz1 + k] = a[i0*ny*nz + j0*nz + k0];
      }
    }
  }
  return x;
}
