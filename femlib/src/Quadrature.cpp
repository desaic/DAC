#include "Quadrature.hpp"

const float Gauss2Pt[2]=
{
 -1.0f/std::sqrt(3.0f),
  1.0f/std::sqrt(3.0f)
};

const float Gauss4Pt[4] =
{
	-std::sqrt( 3.0f / 7.0f + 2.0f * std::sqrt(1.2f) / 7.0f),
	-std::sqrt( 3.0f / 7.0f - 2.0f * std::sqrt(1.2f) / 7.0f),
	 std::sqrt( 3.0f / 7.0f - 2.0f * std::sqrt(1.2f) / 7.0f),
	 std::sqrt( 3.0f / 7.0f + 2.0f * std::sqrt(1.2f) / 7.0f)
};

const float Gauss4Wt[4] =
{
  8*(18.0f - std::sqrt(30.0f)) / 72.0f,
  8*(18.0f + std::sqrt(30.0f)) / 72.0f,
  8*(18.0f + std::sqrt(30.0f)) / 72.0f,
  8*(18.0f - std::sqrt(30.0f)) / 72.0f
};

Quadrature makeGauss2_2D();
Quadrature makeGauss4_2D();
Quadrature makeGauss2();
Quadrature makeNinePt();
Quadrature makeUniform4();

const Quadrature Quadrature::Gauss2_2D = makeGauss2_2D();
const Quadrature Quadrature::Gauss4_2D = makeGauss4_2D();

const Quadrature Quadrature::Gauss2=makeGauss2();
const Quadrature Quadrature::NinePt=makeNinePt();
const Quadrature Quadrature::Uniform4=makeUniform4();

Quadrature makeGauss2_2D()
{
	Quadrature q;
	q.x.resize(4);
  q.w.resize(4, 1.0f);
  q.x[0] = Eigen::Vector3d(Gauss2Pt[0], Gauss2Pt[0], 0);
  q.x[1] = Eigen::Vector3d(Gauss2Pt[1], Gauss2Pt[0], 0);
  q.x[2] = Eigen::Vector3d(Gauss2Pt[0], Gauss2Pt[1], 0);
  q.x[3] = Eigen::Vector3d(Gauss2Pt[1], Gauss2Pt[1], 0);
	return q;
}


Quadrature makeGauss4_2D()
{
  Quadrature q;
  q.x.resize(16);
  q.w.resize(16, 0.25f);
  for (int ii = 0; ii < 4; ii++){
    for (int jj = 0; jj < 4; jj++){
      q.x[ii * 4 + jj] = Eigen::Vector3d(Gauss4Pt[ii], Gauss4Pt[jj],0);
      q.w[ii * 4 + jj] = Gauss4Wt[ii] * Gauss4Wt[jj];
    }
  }
  return q;
}

Quadrature makeGauss2()
{
  Quadrature q;
  q.x.resize(8);
  q.w.resize(8, 1.0f);
  q.x[0] = Eigen::Vector3d(Gauss2Pt[0],Gauss2Pt[0],Gauss2Pt[0]);
  q.x[1] = Eigen::Vector3d(Gauss2Pt[0],Gauss2Pt[0],Gauss2Pt[1]);
  q.x[2] = Eigen::Vector3d(Gauss2Pt[0],Gauss2Pt[1],Gauss2Pt[0]);
  q.x[3] = Eigen::Vector3d(Gauss2Pt[0],Gauss2Pt[1],Gauss2Pt[1]);
  q.x[4] = Eigen::Vector3d(Gauss2Pt[1],Gauss2Pt[0],Gauss2Pt[0]);
  q.x[5] = Eigen::Vector3d(Gauss2Pt[1],Gauss2Pt[0],Gauss2Pt[1]);
  q.x[6] = Eigen::Vector3d(Gauss2Pt[1],Gauss2Pt[1],Gauss2Pt[0]);
  q.x[7] = Eigen::Vector3d(Gauss2Pt[1],Gauss2Pt[1],Gauss2Pt[1]);
  return q;
}

Quadrature makeNinePt()
{
  Quadrature q;
  q.x.resize(9);
  q.w.resize(9, 8/9.0f);
  q.x[0] = Eigen::Vector3d(Gauss2Pt[0], Gauss2Pt[0], Gauss2Pt[0]);
  q.x[1] = Eigen::Vector3d(Gauss2Pt[0], Gauss2Pt[0], Gauss2Pt[1]);
  q.x[2] = Eigen::Vector3d(Gauss2Pt[0], Gauss2Pt[1], Gauss2Pt[0]);
  q.x[3] = Eigen::Vector3d(Gauss2Pt[0], Gauss2Pt[1], Gauss2Pt[1]);
  q.x[4] = Eigen::Vector3d(Gauss2Pt[1], Gauss2Pt[0], Gauss2Pt[0]);
  q.x[5] = Eigen::Vector3d(Gauss2Pt[1], Gauss2Pt[0], Gauss2Pt[1]);
  q.x[6] = Eigen::Vector3d(Gauss2Pt[1], Gauss2Pt[1], Gauss2Pt[0]);
  q.x[7] = Eigen::Vector3d(Gauss2Pt[1], Gauss2Pt[1], Gauss2Pt[1]);
  q.x[8] = Eigen::Vector3d(0,0,0);
  return q;
}

Quadrature makeUniform4()
{
  Quadrature q;
  q.w.resize(64, 8.0f/64.0f);
  int nPt = 4;
  float dx = 2.0f/nPt;
  float x0 = -1;
  for(int ii = 0; ii<nPt; ii++){
    for(int jj = 0; jj<nPt; jj++){
      for(int kk = 0; kk<nPt; kk++){
        Eigen::Vector3d p(x0+dx*(ii+0.5f),x0+dx*(jj+0.5f),x0+dx*(kk+0.5f));
        q.x.push_back(p);
      }
    }
  }
  return q;
}

Quadrature::Quadrature()
{}
