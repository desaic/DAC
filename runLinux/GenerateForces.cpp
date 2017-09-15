#include "GenerateForces.hpp"
#include "FileUtil.hpp"
#include <iostream>
#include <fstream>
#include <vector>
///@brief edges are ordered in increasing order in axis
const int EDGE_AXIS = 4;
const int N_VERT = 8;
int Edges[3][4][2] = {
    { { 0, 4 }, { 1, 5 }, { 2, 6 }, { 3, 7 } },
    { { 0, 2 }, { 1, 3 }, { 4, 6 }, { 5, 7 } },
    { { 0, 1 }, { 2, 3 }, { 4, 5 }, { 6, 7 } }
};

const float compressScale = 1;
int Compress[3][2][4] = {
    { { 0, 1, 2, 3 }, { 4, 5, 6, 7 } },
    { { 0, 1, 4, 5 }, { 2, 3, 6, 7 } },
    { { 0, 2, 4, 6 }, { 1, 3, 5, 7 } }
};

///@brief block has less resistance to twisting.
const float twistScale = 0.5;
int Twist[3][2][4] = {
    { { 0, 4, 6, 2 }, { 1, 5, 7, 3 } },
    { { 0, 1, 5, 4 }, { 3, 7, 6, 2 } },
    { { 0, 1, 3, 2 }, { 4, 5, 7, 6 } }
};

int Shear[3][2][4] = {
    { { 0, 4, 3, 7 }, { 1, 5, 2, 6 } },
    { { 0, 2, 5, 7 }, { 1, 3, 4, 6 } },
    { { 0, 1, 6, 7 }, { 2, 3, 4, 5 } }
};

const float shearForceScale = 0.5;
int ShearForce[3][2][3] =
{
  { { 0, -1, -1 }, { 0, 1, -1 } },
  { { -1, 0, -1 }, { 1, 0, -1 } },
  { { -1, -1, 0 }, { 1, -1, 0 } },
};

const float BendScale = 0.8;
int Bend[2][2][2] =
{
  { { 0, 1 }, { 2, 3 } },
  { { 0, 2 }, { 1, 3 } }
};

const float SurfBendScale = 0.8;
int SurfBend[2][2] =
{
  { 0, 3 }, { 1, 2 }
};

const float RandScale = 1;

int VertInt[8][3] =
{
  { 0, 0, 0 },
  { 0, 0, 1 },
  { 0, 1, 0 },
  { 0, 1, 1 },
  { 1, 0, 0 },
  { 1, 0, 1 },
  { 1, 1, 0 },
  { 1, 1, 1 }
};


//trilinear fit
float arr_efs[27][3] = { { -1, -1, -1 }, { -1, -1, 0 }, { -1, -1, 1 },
{ -1, 0, -1 }, { -1, 0, 0 }, { -1, 0, 1 },
{ -1, 1, -1 }, { -1, 1, 0 }, { -1, 1, 1 },
{ 0, -1, -1 }, { 0, -1, 0 }, { 0, -1, 1 },
{ 0, 0, -1 }, { 0, 0, 0 }, { 0, 0, 1 },
{ 0, 1, -1 }, { 0, 1, 0 }, { 0, 1, 1 },
{ 1, -1, -1 }, { 1, -1, 0 }, { 1, -1, 1 },
{ 1, 0, -1 }, { 1, 0, 0 }, { 1, 0, 1 },
{ 1, 1, -1 }, { 1, 1, 0 }, { 1, 1, 1 } };


void scaleArr(std::vector<float> & aa, float cc)
{
  for (unsigned int ii = 0; ii<aa.size(); ii++){
    aa[ii] *= cc;
  }
}

///@brief normalize L1 sum to 1.
void normalizeSum(std::vector<float> & aa)
{
  float sum = 0;
  for (unsigned int ii = 0; ii<aa.size(); ii++){
    sum += std::abs(aa[ii]);
  }
  scaleArr(aa, 1.0f / (0.00001 + sum));
}

void ForceSample::addBases()
{
  ff = baseForce;
  int idx = baseForce.size();
  for (unsigned int ii = 0; ii<baseForce.size(); ii++){
    ff.push_back(baseForce[ii]);
    for (unsigned int jj = 0; jj<ff[idx].size(); jj++){
      ff[idx][jj] = -ff[idx][jj];
    }
    idx++;
  }
}

///@brief scale forces so that the magnitude of the maximum force is 1
void scaleForce(std::vector<Eigen::Vector3f> & ff)
{
  float maxLen = 0;
  for (unsigned int ii = 0; ii<ff.size(); ii++){
    float len2 = ff[ii].norm();
    if (len2>maxLen){
      maxLen = len2;
    }
  }
  if (maxLen<0.0000001){
    return;
  }
  for (unsigned int ii = 0; ii<ff.size(); ii++){
    ff[ii] /= (maxLen / RandScale);
  }
}

void
ForceSample::saveForces(const char * filename)
{
  FileUtilOut out;
  out.open(filename);

  out.out << ff.size() << "\n" << ff[0].size() << "\n";
  out.out << nMag << "\n";
  for (unsigned int i = 0; i<ff.size(); i++){
    for (unsigned int j = 0; j < ff[i].size(); j++){
      out.out << ff[i][j]<<"\n";
    }
  }

  out.close();
}

void
ForceSample::loadForces(const char * filename)
{
  FileUtilIn in;
  in.open(filename);
  int size;
  in.in >> size;
  ff.resize(size);
  in.in >> size;
  in.in >> nMag;
  for (unsigned int ii = 0; ii<ff.size(); ii++){
    ff[ii].resize(size);
    for (int j = 0; j < size; j++){
      in.in >> ff[ii][j][0] >> ff[ii][j][1] >> ff[ii][j][2];
    }
  }

  in.close();
}

ForceSample::ForceSample() :forceScale(10), nMag(5)
{
  //  generateBase();
}

void ForceSample::generateBase()
{
  generateModeForce();
  generateAxisForce();
}

void ForceSample::scaleForces()
{
  for (size_t i = 0; i < baseForce.size(); i++){
    
    for (int j = 0; j < nMag; j++){
      float scale = (j + 1.0) / nMag * forceScale;
      std::vector<Eigen::Vector3d> fvec(baseForce[i].size());
      for (size_t k = 0; k < baseForce[i].size(); k++){
        fvec[k] = scale * baseForce[i][k];
      }
      ff.push_back(fvec);
    }
    for (int j = 0; j < nMag; j++){
      float scale = (nMag - j) / (float)nMag * forceScale;
      std::vector<Eigen::Vector3d> fvec(baseForce[i].size());
      for (size_t k = 0; k < baseForce[i].size(); k++){
        fvec[k] = scale * baseForce[i][k];
      }
      ff.push_back(fvec);
    }
    for (int j = 0; j < nMag; j++){
      float scale = (j + 1.0) / nMag * forceScale;
      std::vector<Eigen::Vector3d> fvec(baseForce[i].size());
      scale = -scale;
      for (size_t k = 0; k < baseForce[i].size(); k++){
        fvec[k] = scale * baseForce[i][k];
      }
      ff.push_back(fvec);
    }
    for (int j = 0; j < nMag; j++){
      float scale = (nMag-j) / (float)nMag * forceScale;
      std::vector<Eigen::Vector3d> fvec(baseForce[i].size());
      scale = -scale;
      for (size_t k = 0; k < baseForce[i].size(); k++){
        fvec[k] = scale * baseForce[i][k];
      }
      ff.push_back(fvec);
    }
  }
}

void ForceSample::generateAxisForce()
{
  for (int axis = 0; axis<3; axis++){
    for (int jj = 0; jj<EDGE_AXIS; jj++){
      std::vector<Eigen::Vector3d> forces;
      forces.resize(N_VERT, Eigen::Vector3d::Zero());
      forces[Edges[axis][jj][0]][axis] = 1;
      forces[Edges[axis][jj][1]][axis] = -1;
      baseForce.push_back(forces);
    }
  }
}

void ForceSample::generateTwist(bool flip)
{
  for (int axis = 0; axis<3; axis++){
    std::vector<Eigen::Vector3d> forces;
    forces.resize(N_VERT, Eigen::Vector3d::Zero());
    for (int ii = 0; ii<4; ii++){
      int vi = Twist[axis][0][ii];
      int next, prev;
      //      if(flip){
      //        next =Twist[axis][0][(ii+1)%4];
      //        prev = Twist[axis][1][(ii+3)%4];
      //      }else{
      next = Twist[axis][0][(ii + 3) % 4];
      prev = Twist[axis][0][(ii + 1) % 4];
      //      }
      for (int jj = 0; jj<3; jj++){
        forces[vi][jj] = VertInt[next][jj] - VertInt[prev][jj];
      }

      vi = Twist[axis][1][ii];
      next = Twist[axis][1][(ii + 3) % 4];
      prev = Twist[axis][1][(ii + 1) % 4];

      for (int jj = 0; jj<3; jj++){
        forces[vi][jj] = VertInt[prev][jj] - VertInt[next][jj];
      }
    }

    for (unsigned int ii = 0; ii<forces.size(); ii++){
      forces[ii] *= twistScale;
    }
    baseForce.push_back(forces);
  }
}

void ForceSample::generateModeForce()
{
  //compression
  for (int axis = 0; axis<3; axis++){
    std::vector<Eigen::Vector3d> forces;
    forces.resize(N_VERT, Eigen::Vector3d::Zero());
    for (int ii = 0; ii<4; ii++){
      forces[Compress[axis][0][ii]][axis] = compressScale;
      forces[Compress[axis][1][ii]][axis] = -compressScale;
    }
    baseForce.push_back(forces);
  }
  //shear
  for (int axis = 0; axis<3; axis++){
    std::vector<Eigen::Vector3d> forces;
    forces.resize(N_VERT, Eigen::Vector3d::Zero());
    for (int kk = 0; kk<2; kk++){
      for (int ll = 0; ll<3; ll++){
        forces[Shear[axis][0][kk]][ll] = shearForceScale * ShearForce[axis][0][ll];
        forces[Shear[axis][0][kk + 2]][ll] = -shearForceScale * ShearForce[axis][0][ll];
        forces[Shear[axis][1][kk]][ll] = shearForceScale * ShearForce[axis][1][ll];
        forces[Shear[axis][1][kk + 2]][ll] = -shearForceScale * ShearForce[axis][1][ll];
      }
    }
    baseForce.push_back(forces);
  }

  //bending
  for (int axis = 0; axis<3; axis++){
    for (int jj = 0; jj<2; jj++){
      std::vector<Eigen::Vector3d> forces;
      forces.resize(N_VERT, Eigen::Vector3d::Zero());

      for (int kk = 0; kk<2; kk++){
        int eIdx = Bend[jj][0][kk];
        int v0 = Edges[axis][eIdx][0];
        int v1 = Edges[axis][eIdx][1];
        for (int ll = 0; ll<3; ll++){
          forces[v0][ll] = VertInt[v0][ll] - VertInt[v1][ll];
          forces[v1][ll] = -forces[v0][ll];
        }

        eIdx = Bend[jj][1][kk];
        v0 = Edges[axis][eIdx][0];
        v1 = Edges[axis][eIdx][1];
        for (int ll = 0; ll<3; ll++){
          forces[v0][ll] = VertInt[v1][ll] - VertInt[v0][ll];
          forces[v1][ll] = -forces[v0][ll];
        }
      }
      for (auto kk = 0; kk<forces.size(); kk++){
        forces[kk] *= BendScale;
      }
      baseForce.push_back(forces);
    }
  }
  generateTwist();
}
