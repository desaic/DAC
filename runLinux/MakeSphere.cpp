#include <Eigen/Dense>
#include "ArrayUtil.hpp"
#include "ElementMesh.hpp"
#include "FileUtil.hpp"
struct Quad
{
  ///@brief center of quad in spherical coord
  ///theta is angle from y axis, phi is on x-z plane
  ///theta in [0, PI/2], phi in [0, 2PI)
  Eigen::Vector3d s_center;
  //maximum 4 neighbors
  int nbr[4];
  int v[4];
  Quad(){
    for (int i = 0; i < 4; i++){
      nbr[i] = -1;
      v[i] = -1;
    }
  }

  void addNbr(int n){
    bool dup = false;
    int last = -1;
    for (int i = 0; i < 4; i++){
      if (nbr[i] == n){
        dup = true;
        break;
      }
      if (nbr[i] < 0){
        last = i;
        break;
      }
    }
    if (!dup){
      if (last >= 0){
        nbr[last] = n;
      }
      else{
        std::cout << "Add nbr too many nbrs.\n";
      }
    }
  }

};

struct SphereMesh
{
  float r;
  float circleDist;
  float cubeSize;
  int tipIdx;
  std::vector<Eigen::VectorXi> cubes;
  std::vector<Eigen::Vector3d> cubev;
  std::vector<Quad> quads;
  std::vector<Eigen::Vector3d> v;
  //indices of quads in the cross
  std::vector<int> cross[4];
  ///@brief list of circles. Each circle has 4 arcs. Each arc has a list of indices in increasing
  ///phi order.
  std::vector<std::vector<std::vector<int> > >circle;

  std::vector<std::vector<std::vector<int> > >diags;

  SphereMesh() :r(0), tipIdx(-1){}
};

///@param dir 0 for vertical theta and 1 for horizontal phi
///@param idx indices of the edge including end points in order.
///@param angle angle of the other param on sphere.
///@brief q1 should have larger or equal phi to avoid wrap around ambiguity.
void connectQuads(SphereMesh & m, int q0idx, int q1idx, int nCube, int dir, float angle,
  std::vector<int> & idx)
{
  idx.push_back(q0idx);
  if (dir == 0){
    for (int i = 0; i < nCube; i++){
      Quad & q0 = m.quads[q0idx];
      Quad & q1 = m.quads[q1idx];
      double phi = angle;
      double t0 = q0.s_center[1];
      double t1 = q1.s_center[1];

      float theta = t0 + (i + 1) * (t1 - t0) / (nCube + 1);
      Quad qi;
      int qidx = m.quads.size();
      qi.s_center = Eigen::Vector3d(m.r, theta, phi);
      if (i == 0){
        qi.addNbr(q0idx);
        q0.addNbr(qidx);
      }
      else{
        qi.addNbr(qidx - 1);
        m.quads[qidx - 1].addNbr(qidx);
      }
      if (i == nCube - 1){
        qi.addNbr(q1idx);
        q1.addNbr(qidx);
      }
      if (idx.size()>20){
        std::cout << idx[20] << "\n";
      }
      m.quads.push_back(qi);
      idx.push_back(qidx);
    }
  }
  else{
    for (int i = 0; i < nCube; i++){
      Quad & q0 = m.quads[q0idx];
      Quad & q1 = m.quads[q1idx];
      double p0 = q0.s_center[2];
      double p1 = q1.s_center[2];
      if (p1 < p0){
        p1 = p1 + 2 * M_PI;
      }
      double t0 = angle;

      float phi = p0 + (i + 1) * (p1 - p0) / (nCube + 1);
      Quad qi;
      int qidx = m.quads.size();
      qi.s_center = Eigen::Vector3d(m.r, t0, phi);
      if (i == 0){
        qi.addNbr(q0idx);
        q0.addNbr(qidx);
      }
      else{
        qi.addNbr(qidx - 1);
        m.quads[qidx - 1].addNbr(qidx);
      }
      if (i == nCube - 1){
        qi.addNbr(q1idx);
        q1.addNbr(qidx);
      }
      m.quads.push_back(qi);
      idx.push_back(qidx);
    }
  }
  idx.push_back(q1idx);
}

double arclenTheta(double r, double t0, double t1)
{
  return r * std::abs(t1 - t0);
}

double arclenPhi(double r, double theta, double phi0, double phi1)
{
  return r * std::abs(std::sin(theta) * (phi1 - phi0));
}

Eigen::Vector3d sphere2Cart(const Eigen::Vector3d & sp)
{
  Eigen::Vector3d v;
  v[0] = sp[0] * std::sin(sp[1]) * std::cos(sp[2]);
  v[1] = sp[0] * std::cos(sp[1]);
  v[2] = sp[0] * std::sin(sp[1]) * std::sin(sp[2]);
  return v;
}

void writeQuadCenter(const SphereMesh & m, std::string outname)
{
  FileUtilOut out(outname);
  out.out << m.quads.size() << "\n";
  for (size_t i = 0; i < m.quads.size(); i++){
    Eigen::Vector3d x = sphere2Cart(m.quads[i].s_center);
    out.out << x[0] << " " << x[1] << " " << x[2] << "\n";
  }
  out.close();
}

void makeCircle(SphereMesh & m, double dist)
{
  int q0idx, q1idx;
  int crossQuadIdx = -1;
  for (int i = 0; i < (int)m.cross[0].size(); i++){
    double theta = m.quads[m.cross[0][i]].s_center[1];
    if (arclenTheta(m.r, theta, m.quads[m.tipIdx].s_center[1]) > dist){
      crossQuadIdx = i;
      break;
    }
  }
  if (crossQuadIdx < 0){
    std::cout << "Warning using last cube on cross.\n";
    crossQuadIdx = (int)m.cross[0].size() - 1;
  }

  //each circle has 4 arcs
  double theta = m.quads[m.cross[0][crossQuadIdx]].s_center[1];
  q0idx = m.cross[0][crossQuadIdx];
  double phi0 = m.quads[q0idx].s_center[2];
  q1idx = m.cross[2][crossQuadIdx];
  double phi1 = m.quads[q1idx].s_center[2];
  int nCube = arclenPhi(m.r, theta, phi0, phi1) / m.cubeSize;
  //need odd number of cubes for diagonal spokes.
  if (nCube % 2 == 0){
    nCube++;
  }
  std::vector<std::vector<int> > circle(4);
  connectQuads(m, q0idx, q1idx, nCube, 1, theta, circle[0]);
  q0idx = m.cross[2][crossQuadIdx];
  q1idx = m.cross[1][crossQuadIdx];
  connectQuads(m, q0idx, q1idx, nCube, 1, theta, circle[1]);
  q0idx = m.cross[1][crossQuadIdx];
  q1idx = m.cross[3][crossQuadIdx];
  connectQuads(m, q0idx, q1idx, nCube, 1, theta, circle[2]);
  q0idx = m.cross[3][crossQuadIdx];
  q1idx = m.cross[0][crossQuadIdx];
  connectQuads(m, q0idx, q1idx, nCube, 1, theta, circle[3]);
  m.circle.push_back(circle);
}

///@param circIdx bottom circle index
void makeDiag(SphereMesh & m, int circIdx)
{
  //circle has 4 arcs.
  std::vector<std::vector< int > > diags(4);
  for (int i = 0; i < 4; i++){
    int quadInCirc = (int)m.circle[circIdx - 1][i].size() / 2;
    int q0idx = m.circle[circIdx - 1][i][quadInCirc];
    quadInCirc = (int)m.circle[circIdx][i].size() / 2;
    int q1idx = m.circle[circIdx][i][quadInCirc];
    double phi = m.quads[q0idx].s_center[2];
    double t0 = m.quads[q0idx].s_center[1];
    double t1 = m.quads[q1idx].s_center[1];
    int nCube = (int)(arclenTheta(m.r, t0, t1) / m.cubeSize) + 1;
    connectQuads(m, q0idx, q1idx, nCube, 0, phi, diags[i]);
    m.diags.push_back(diags);
  }
}

void makeQuads(SphereMesh & m)
{
  //counter-clockwise
  double localVert[4][2] = { { 0.5, 0.5 }, { -0.5, 0.5 }, { -0.5, -0.5 }, { 0.5, -0.5 } };
  for (size_t i = 0; i < m.quads.size(); i++){
    for (int vi = 0; vi < 4; vi++){

      //sphere center to quad center
      Eigen::Vector3d rvec = sphere2Cart(m.quads[i].s_center);
      Eigen::Vector3d up(0, 1, 0);
      //x axis align with horizontal
      Eigen::Vector3d xa = up.cross(rvec);
      //local quad y axis aligns with vertical
      Eigen::Vector3d ya;
      double eps = 1e-20;
      //special top point. x aligns with x x-axis in world
      if (xa.squaredNorm() < eps || m.quads[i].s_center[1] < eps){
        xa = Eigen::Vector3d(1, 0, 0);
        ya = Eigen::Vector3d(0, 0, -1);
      }
      else{
        xa.normalize();
        ya = rvec.cross(xa).normalized();
      }
      //vertex pos if not exist.
      Eigen::Vector3d x = rvec + m.cubeSize * (localVert[vi][0] * xa + localVert[vi][1] * ya);
      bool hasNbr = false;
      for (int ni = 0; ni < 4; ni++){
        int nbrIdx = m.quads[i].nbr[ni];
        if (nbrIdx < 0){
          break;
        }
        if (m.quads[nbrIdx].v[0] < 0){
          continue;
        }
        //there is a vertex in neighbor for this vertex.
        ///vertex to quad center
        Eigen::Vector3d v0 = rvec - x;
        Eigen::Vector3d nbrCenter = sphere2Cart(m.quads[nbrIdx].s_center);
        Eigen::Vector3d v1 = nbrCenter - x;
        Eigen::Vector3d centerDiff = nbrCenter - rvec;
        //if the two vectors are opposite, the vertex is shared in the middle.
        if (v0.dot(centerDiff) * v1.dot(centerDiff) < 0){
          hasNbr = true;
          //find nearest vertex in nbr
          double minDist = m.r * m.r;
          int minIdx = -1;
          for (int nv = 0; nv < 4; nv++){
            int nvIdx = m.quads[nbrIdx].v[nv];
            Eigen::Vector3d nbrx = m.v[nvIdx];
            double dist = (nbrx - x).squaredNorm();
            if (dist < minDist){
              minDist = dist;
              minIdx = nvIdx;
            }
          }
          m.quads[i].v[vi] = minIdx;
          break;
        }
      }
      if (!hasNbr){
        m.quads[i].v[vi] = m.v.size();
        m.v.push_back(x);
      }
    }
  }
}

void saveQuads(const SphereMesh & m, std::string filename)
{
  FileUtilOut out(filename);
  for (size_t i = 0; i < m.v.size(); i++){
    out.out << "v " << m.v[i][0] << " " << m.v[i][1] << " " << m.v[i][2] << "\n";
  }
  for (size_t i = 0; i < m.quads.size(); i++){
    out.out << "f";
    for (int j = 0; j < 4; j++){
      out.out << " " << (m.quads[i].v[j] + 1);
    }
    out.out << "\n";
  }
  out.close();
}

void makeCubes(SphereMesh & m)
{
  m.cubev.resize(2 * m.v.size());
  for (int i = 0; i < m.v.size(); i++){
    m.cubev[2 * i] = m.v[i];
    m.cubev[2 * i + 1] = m.v[i] - m.cubeSize * m.v[i].normalized();
  }
  m.cubes.resize(m.quads.size(), Eigen::VectorXi(8));
  for (int i = 0; i < (int)m.cubes.size(); i++){
    m.cubes[i][0] = 2 * m.quads[i].v[0] + 1;
    m.cubes[i][1] = 2 * m.quads[i].v[0];
    m.cubes[i][2] = 2 * m.quads[i].v[3] + 1;
    m.cubes[i][3] = 2 * m.quads[i].v[3];
    m.cubes[i][4] = 2 * m.quads[i].v[1] + 1;
    m.cubes[i][5] = 2 * m.quads[i].v[1];
    m.cubes[i][6] = 2 * m.quads[i].v[2] + 1;
    m.cubes[i][7] = 2 * m.quads[i].v[2];
  }
}

void saveCubes(const SphereMesh & m, std::string filename)
{
  //scale from cm to m.
  float scale = 0.01;
  FileUtilOut out(filename);
  out.out << "#verts " << m.cubev.size() << "\n";
  out.out << "#elts " << m.cubes.size() << "\n";
  for (size_t i = 0; i < m.cubev.size(); i++){
    Eigen::Vector3d v = m.cubev[i];
    v *= scale;
    out.out << v[0] << " " << v[1] << " " << v[2] << "\n";
  }
  for (size_t i = 0; i < m.cubes.size(); i++){
    out.out << "8";
    for (int j = 0; j < 8; j++){
      out.out << " " << (m.cubes[i][j]);
    }
    out.out << "\n";
  }
  out.out << "parts\n" << m.cubes.size() << "\n";
  for (size_t i = 0; i < m.cubes.size(); i++){
    out.out << "1\n";
  }
  out.out << "depth\n" << m.cubes.size() << "\n";
  for (size_t i = 0; i < m.cubes.size(); i++){
    out.out << "1\n";
  }

  out.close();
}

void makeSphere()
{
  //cm
  float r = 9;
  float cubeSize = 0.25;
  SphereMesh m;
  int nCirc = 4;
  m.circleDist = r * 0.5 * M_PI / nCirc;
  m.cubeSize = cubeSize;
  m.r = r;
  //make cross
  float arcLen = 0.5 * r * M_PI;
  int nCube = (int)(arcLen / cubeSize);
  float phi = M_PI;
  //along x axis
  //make endpoint quads
  Quad q0, q1, qtop;
  qtop.s_center = Eigen::Vector3d(r, 0, 0);
  q0.s_center = Eigen::Vector3d(r, M_PI / 2, M_PI);
  q1.s_center = Eigen::Vector3d(r, M_PI / 2, 0);
  m.quads.push_back(qtop);
  m.tipIdx = 0;
  m.quads.push_back(q0);
  m.quads.push_back(q1);
  connectQuads(m, m.tipIdx, 1, nCube, 0, q0.s_center[2], m.cross[0]);
  connectQuads(m, m.tipIdx, 2, nCube, 0, q1.s_center[2], m.cross[1]);

  //along z axis
  q0.s_center = Eigen::Vector3d(r, M_PI / 2, 1.5 * M_PI);
  q1.s_center = Eigen::Vector3d(r, M_PI / 2, 0.5 * M_PI);
  int q0idx = m.quads.size();
  m.quads.push_back(q0);
  int q1idx = q0idx + 1;
  m.quads.push_back(q1);
  connectQuads(m, m.tipIdx, q0idx, nCube, 0, q0.s_center[2], m.cross[2]);
  connectQuads(m, m.tipIdx, q1idx, nCube, 0, q1.s_center[2], m.cross[3]);

  //make first circle with odd number of cubes on each quater arc
  //find the quad on the cross circleDist arc distance away from top
  makeCircle(m, m.circleDist);
  makeCircle(m, 2 * m.circleDist);
  makeCircle(m, 3 * m.circleDist);
  makeCircle(m, 4 * m.circleDist);

  //make diagonal spokes.
  for (int i = 1; i < (int)m.circle.size(); i++){
    makeDiag(m, i);
  }
  //writeQuadCenter(m, "qc.txt");
  makeQuads(m);
  //saveQuads(m, "sphereQuads.obj");
  makeCubes(m);
  saveCubes(m, "cubes.txt");
}

void deformToSphere(ElementMesh * em)
{
  //place at center of mass for x y, 0 for z
  Eigen::Vector3d centerOfMass = em->centerOfMass();
  Eigen::Vector3d mn, mx;
  BBox(em->X, mn, mx);
  for (size_t i = 0; i < em->X.size(); i++){
    em->X[i][0] -= centerOfMass[0];
    em->X[i][1] -= centerOfMass[1];
    em->X[i][2] -= mn[2];
  }

  float d0 = 0.095;
  float d1 = 0.1;
  float r0 = d0 / 2;
  float r1 = d1 / 2;
  em->x = em->X;

  float thickness = mx[2] - mn[2];
  for (size_t i = 0; i < em->X.size(); i++){
    float r0w = (mx[2] - em->X[i][2]) / thickness;
    float r = r0 * r0w + r1*(1 - r0w);
    double phi = std::atan2(em->X[i][1], em->X[i][0]);
    double xydist = em->X[i][0] * em->X[i][0] + em->X[i][1] * em->X[i][1];
    xydist = std::sqrt(xydist);
    if (xydist > r1){
      xydist = r1;
    }
    double theta = std::asin(xydist / r1);

    em->X[i][0] = r * std::sin(theta) * std::cos(phi);
    em->X[i][1] = r * std::sin(theta) * std::sin(phi);
    em->X[i][2] = r * std::cos(theta);
  }
  em->x = em->X;
}
