#include <iostream>
#include <set>

#include "ArrayUtil.hpp"
#include "EigenUtil.hpp"
#include "Element.hpp"
#include "ElementHex.hpp"
#include "ElementMeshUtil.hpp"
#include "ElementMesh.hpp"
#include "FileUtil.hpp"
#include "MaterialQuad.hpp"
#include "TrigMesh.hpp"
#include "Quadrature.hpp"

typedef std::map<int, int> IntMap;

int HexFaces[6][4] = {
    { 0, 1, 3, 2 }, { 4, 6, 7, 5 },
    { 0, 4, 5, 1 }, { 2, 3, 7, 6 },
    { 0, 2, 6, 4 }, { 1, 5, 7, 3 } };
static const int sw[8][3] =
{ { -1, -1, -1 },
{ -1, -1, 1 },
{ -1, 1, -1 },
{ -1, 1, 1 },
{ 1, -1, -1 },
{ 1, -1, 1 },
{ 1, 1, -1 },
{ 1, 1, 1 }
};
struct IntGrid{

  //what indices in each cell.
  std::vector < std::vector<int> > g;

  //points added to the grid.
  std::vector<Eigen::Vector3d> pts;
  
  //grid size.
  Eigen::Vector3i s;
  
  Eigen::Vector3d origin;
  
  //grid size.
  Eigen::Vector3d dx;

  void allocate(Eigen::Vector3i size){
    s = size;
    int nCell = size[0] * size[1] * size[2];
    g.resize(nCell);
  }

  int grid2Linear(const Eigen::Vector3i & i){
    return i[0] * s[1] * s[2] + i[1] * s[2] + i[2];
  }

  Eigen::Vector3i linear2Grid(int i){
    std::vector<int> idx(3);
    idx[0] = i / (s[2] * s[1]);
    i -= idx[0] * (s[2] * s[1]);
    idx[1] = i / s[2];
    i -= idx[1] * s[2];
    idx[2] = i;
  }

  Eigen::Vector3i find(const Eigen::Vector3d & x){
    Eigen::Vector3d rel = x - origin;
    Eigen::Vector3i idx;
    for (int i = 0; i < 3; i++){
      int j = (int)(rel[i] / dx[i]);
      if (j < 0){
        j = 0;
      }
      if (j >= s[i]){
        j = s[i] - 1;
      }
      idx[i] = j;
    }
    return idx;
  }

  void addPoint(const Eigen::Vector3d & x){
    Eigen::Vector3i idx = find(x);
    int ptIdx = pts.size();
    pts.push_back(x);
    int gIdx = grid2Linear(idx);
    g[gIdx].push_back(ptIdx);
  }

  bool inBound(const Eigen::Vector3i  & idx){
    for (int i = 0; i < 3; i++){
      if (idx[i] < 0 || idx[i] >= s[i]){
        return false;
      }
    }
    return true;
  }

  ///@return -1 if no point near 1-ring.
  int findNbr(const Eigen::Vector3d & x){
    
    Eigen::Vector3i centerIdx = find(x);
    int minIdx = -1;
    double minDist = 0;
    for (int i = -1; i < 2; i++){
      for (int j = -1; j < 2; j++){
        for (int k = -1; k < 2; k++){
          Eigen::Vector3i idx;
          idx[0] = centerIdx[0] + i;
          idx[1] = centerIdx[1] + j;
          idx[2] = centerIdx[2] + k;
          if (inBound(idx)){
            int nbrIdx = grid2Linear(idx);
            for (size_t l = 0; l < g[nbrIdx].size(); l++){
              int x1idx = g[nbrIdx][l];
              Eigen::Vector3d x1 = pts[x1idx];
              double squaredDist = (x1 - x).squaredNorm();
              if (minIdx < 0 || squaredDist < minDist){
                minDist = squaredDist;
                minIdx = x1idx;
              }
            }
          }
        }
      }
    }
    return minIdx;
  }
};

void equalityConstraint(std::vector<int> dofs, ElementMesh * em,
  Eigen::SparseMatrix<double> & K)
{
  if (dofs.size() < 2){
    return;
  }
  int nConstraint = (int)dofs.size() - 1;
  int N = (int)K.rows();
  int N1 = N + nConstraint;
  float constraintCoef = 1000;
  K.conservativeResize(N1, N1);
  Eigen::SparseMatrix<double> KConstraint(N1, N1);
  //each constraint affects two unknowns
  KConstraint.reserve(Eigen::VectorXi::Constant(N1, 2));
  for (int i = 1; i < dofs.size(); i++){
    int row = N + i - 1;
    KConstraint.insert(row, dofs[0]) = constraintCoef;
    KConstraint.insert(row, dofs[i]) = -constraintCoef;
    //diagonal must exist for pardiso
    KConstraint.insert(row, row) = 0;
  }

  //only lower triangle of the matrix is actually used in static solver
  K += KConstraint;
}

int grid2Linear(const Eigen::Vector3i & i, int N)
{
  return i[0]*N*N + i[1]*N + i[2];
}

Eigen::Vector3d trilinear(Eigen::Vector3d p, std::vector<Eigen::Vector3d> & x)
{
  Eigen::Vector3d v = Eigen::Vector3d::Zero();
  for (size_t i = 0; i < x.size(); i++){
    v += (1.0f / 8) * (1 + sw[i][0] * p[0])
      *(1 + sw[i][1] * p[1]) *(1 + sw[i][2] * p[2]) * x[i];
  }
  return v;
}

ElementMesh * subdivideHex(ElementMesh * em, int nSub)
{
  
  int nFineV = (nSub + 1) * (nSub + 1) * (nSub + 1);
  //natural coordinate of fine vertices.
  std::vector<Eigen::Vector3d> natCoord(nFineV);
  for (int i = 0; i <= nSub; i++){
    for (int j = 0; j <= nSub; j++){
      for (int k = 0; k <= nSub; k++){
        int idx = i * (nSub + 1)* (nSub + 1) + j * (nSub + 1) + k;
        natCoord[idx] = Eigen::Vector3d(i*2.0 / nSub - 1, j*2.0 / nSub - 1, k*2.0 / nSub - 1);
      }
    }
  }
  //3d index of first vertex in a cube.
  int nFineE = nSub * nSub * nSub ;
  std::vector<Eigen::Vector3i> v0idx(nFineE);
  for (int i = 0; i < nSub; i++){
    for (int j = 0; j < nSub; j++){
      for (int k = 0; k < nSub; k++){
        int idx = i * (nSub )* (nSub ) + j * (nSub ) + k;
        v0idx[idx] = Eigen::Vector3i(i, j, k);
      }
    }
  }
  int oneV[8][3] = { { 0, 0, 0 }, { 0, 0, 1 }, { 0, 1, 0 }, { 0, 1, 1 },
  { 1, 0, 0 }, { 1, 0, 1 }, { 1, 1, 0 }, { 1, 1, 1 } };
  //local vertex indices of fine elements.
  std::vector<std::vector<int> > fineVi(nFineE);
  //loop over fine elements
  for (int i = 0; i < nFineE; i++){
    fineVi[i].resize(8);
    //loop over 8 vertices of each fine element.
    for (int j = 0; j < 8; j++){
      Eigen::Vector3i vlocalIdx;
      for (int k = 0; k < 3; k++){
        vlocalIdx[k] = oneV[j][k];
      }
      vlocalIdx += v0idx[i];
      fineVi[i][j] = grid2Linear(vlocalIdx, nSub+1);
    }
  }

  std::vector<Element * > e_new;
  IntGrid grid;
  Eigen::Vector3d mn, mx;
  BBox(em->X, mn, mx);
  grid.origin = mn;
  double cubeSize = (1.0/nSub) * (em->X[0] - em->X[1]).norm();
  double eps = cubeSize * 0.1;
  for (int i = 0; i < 3; i++){
    grid.s[i] = (int)((mx[i] - mn[i]) / cubeSize )+ 1;
    grid.dx[i] = cubeSize;
  }
  grid.allocate(grid.s);

  //make new vertices with duplicates
  for (size_t i = 0; i < em->e.size(); i++){
    std::vector<int> newvidx(natCoord.size(), -1);
    Element * ele = em->e[i];
    std::vector<Eigen::Vector3d> x(8);
    for (int j = 0; j < 8; j++){
      x[j] = em->X[ele->at(j)];
    }
    //loop over fine vertices
    for (int j = 0; j < (int)natCoord.size(); j++){
      Eigen::Vector3d vfine = Eigen::Vector3d::Zero();
      vfine = trilinear(natCoord[j], x);
      int ni = grid.findNbr(vfine);
      int vj = grid.pts.size();
      if (ni >= 0){
        Eigen::Vector3d vnear = grid.pts[ni];
        double dist = (vnear - vfine).norm();
        if (dist > eps){
          grid.addPoint(vfine);
        }
        else{
          vj = ni;
        }
      }
      else{
        grid.addPoint(vfine);
      }

      newvidx[j] = vj;
    }
    //loop over fine elements
    for (int j = 0; j < (int)fineVi.size(); j++){
      ElementHex * ele = new ElementHex();
      for (int k = 0; k < fineVi[j].size(); k++){
        (*ele)[k] = newvidx[fineVi[j][k]];
      }
      e_new.push_back(ele);
    }
  }
  ElementMesh * em_new = new ElementMesh();
  em_new->e = e_new;
  em_new->X = grid.pts;
  return em_new;
  //FileUtilOut out("fine.txt");
  //em_new->initArrays();
  //em_new->saveMesh(out.out);
  //TrigMesh tm;
  //hexToTrigMesh(em_new, &tm);
  //tm.save_obj("fine.obj");
  //out.close();
  //saveAbq(em_new, "sphere2.abq");
  //delete em_new;
}

void specialVertexBound(ElementMesh * m, float lb)
{
  m->topVerts.clear();
  if (m->dim == 2){
    //find top left element
    std::vector<std::vector< int> > grid;
    m->grid2D(grid);
    int nx = (int)grid.size();
    int ny = (int)grid[0].size();
    int topi = -1, topj = -1;
    for (int jj = ny - 1; jj >= 0; jj--){
      for (int ii = 0; ii<nx; ii++){
        if (grid[ii][jj] >= 0){
          topi = ii;
          topj = jj;
          break;
        }
      }
      if (topi>0){
        break;
      }
    }

    if (topi<0){
      std::cout << "Error initializing contact verts. Empty grid\n";
      return;
    }

    int eidx = grid[topi][topj];
    //top left vertex
    m->topVerts.push_back(m->e[eidx]->at(1));
    int ii = topi + 1;
    //find 1 element past right most element
    for (; ii<nx; ii++){
      if (grid[ii][topj]<0){
        break;
      }
    }
    //top right vertex
    eidx = grid[ii - 1][topj];
    m->topVerts.push_back(m->e[eidx]->at(3));

    int jj = topj - 1;
    //find downward neighbor of top left element until break.
    for (; jj >= 0; jj--){
      if (grid[topi][jj]<0){
        break;
      }
    }

    eidx = grid[topi][jj + 1];

    m->specialVert = m->e[eidx]->at(0);
    std::cout << "special vert " << m->specialVert << "\n";
    std::cout << lb << "\n";
    m->lb[m->specialVert][1] = lb;
  }
  else if (m->dim == 3){
    //find top left front element
    std::vector<std::vector< std::vector<int> > > grid;
    m->grid3D(grid);
    int nx = (int)grid.size();
    int ny = (int)grid[0].size();
    int nz = (int)grid[0][0].size();
    int topi = -1, topj = -1, topk = -1;
    for (int jj = ny - 1; jj >= 0; jj--){
      for (int ii = 0; ii<nx; ii++){
        for (int kk = 0; kk<nz; kk++){
          if (grid[ii][jj][kk] >= 0){
            topi = ii;
            topj = jj;
            topk = kk;
            break;
          }
        }
        if (topi >= 0){
          break;
        }
      }
      if (topi >= 0){
        break;
      }
    }

    if (topi<0){
      std::cout << "Error initializing contact verts. Empty grid\n";
      return;
    }

    int eidx = grid[topi][topj][topk];
    //top left vertex
    m->topVerts.push_back(m->e[eidx]->at(2));
    int ii = topi + 1;
    //find 1 element past left most element
    for (; ii<nx; ii++){
      if (grid[ii][topj][topk]<0){
        break;
      }
    }
    //find top right vertex used for measuring top angle
    eidx = grid[ii - 1][topj][topk];
    m->topVerts.push_back(m->e[eidx]->at(6));

    int jj = topj - 1;
    //find downward neighbor of top left element until break.
    for (; jj >= 0; jj--){
      if (grid[topi][jj][topk]<0){
        break;
      }
    }

    eidx = grid[topi][jj + 1][topk];
    m->specialVert = m->e[eidx]->at(0);
    std::cout << "special vert " << m->specialVert << "\n";
    std::cout << lb << "\n";
    for (int kk = topk; kk<nz; kk++){
      eidx = grid[topi][jj + 1][kk];
      if (eidx<0){
        continue;
      }
      m->lb[m->e[eidx]->at(0)][1] = lb;
      m->lb[m->e[eidx]->at(1)][1] = lb;
    }
  }
}

void updateSurfVert(const std::vector<Eigen::Vector3d> & x, TrigMesh * tm, const std::vector<int> & vidx)
{
  for (unsigned int ii = 0; ii<vidx.size(); ii++){
    tm->v[ii] = x[vidx[ii]];
  }
}

//find neighboring elements of a vertex
void elementNeighbors(ElementMesh * em,
  std::vector<std::vector< int> > & eleNeighbor)
{
  eleNeighbor.clear();
  eleNeighbor.resize(em->X.size());
  for (unsigned int ii = 0; ii<em->e.size(); ii++){
    for (int jj = 0; jj<em->e[ii]->nV(); jj++){
      eleNeighbor[em->e[ii]->at(jj)].push_back(ii);
    }
  }
}

int findHexFace(const IntMap & m, int ei, ElementMesh * em){
  for (int fi = 0; fi<6; fi++){
    int cnt = 0;
    for (int fv = 0; fv<4; fv++){
      if (m.at(em->e[ei]->at(HexFaces[fi][fv])) == 2){
        cnt++;
      }
    }
    if (cnt == 4){
      return fi;
    }
  }
  return 0;
}

void hexToTrigMesh(ElementMesh * em, TrigMesh * surf)
{
  std::vector<std::vector<bool> > exterior(em->e.size());
  for (unsigned int ii = 0; ii<exterior.size(); ii++){
    unsigned int nV = em->e[ii]->nV();
    exterior[ii].resize(nV, true);
  }

  //mark exterior faces
  std::vector<std::vector<int > > eleNeighbor;
  elementNeighbors(em, eleNeighbor);
  for (unsigned int ii = 0; ii<em->X.size(); ii++){
    for (unsigned int nj = 0; nj<eleNeighbor[ii].size(); nj++){
      int jj = eleNeighbor[ii][nj];
      for (unsigned int nk = nj + 1; nk<eleNeighbor[ii].size(); nk++){
        int kk = eleNeighbor[ii][nk];
        IntMap im;
        for (int vv = 0; vv<em->e[jj]->nV(); vv++){
          im[em->e[jj]->at(vv)] = 1;
        }
        for (int vv = 0; vv<em->e[kk]->nV(); vv++){
          int key = em->e[kk]->at(vv);
          auto entry = im.find(key);
          if (entry != im.end()){
            im[key] = 2;
          }
          else{
            im[key] = 1;
          }
        }
        if (im.size() == 12){
          //share a face
          int fi = findHexFace(im, jj, em);
          exterior[jj][fi] = false;
          fi = findHexFace(im, kk, em);
          exterior[kk][fi] = false;
        }
      }
    }
  }
  //save exterior vertices and exterior faces
  std::vector<int> e2surf(em->X.size(), -1);
  int vCnt = 0;
  //save vertices
  surf->v.clear();
  surf->t.clear();
  surf->vidx.clear();
  for (unsigned int ii = 0; ii<em->e.size(); ii++){
    for (int fi = 0; fi<6; fi++){
      if (!exterior[ii][fi]){
        continue;
      }
      for (int fv = 0; fv<4; fv++){
        int vi = em->e[ii]->at(HexFaces[fi][fv]);
        if (e2surf[vi]<0){
          e2surf[vi] = vCnt;
          surf->v.push_back(em->x[vi]);
          surf->vidx.push_back(vi);
          vCnt++;
        }
      }
    }
  }
  //save faces
  for (unsigned int ii = 0; ii<em->e.size(); ii++){
    for (int fi = 0; fi<6; fi++){
      if (!exterior[ii][fi]){
        continue;
      }
      int trigs[2][3] = { { 0, 1, 2 }, { 2, 3, 0 } };
      for (int ti = 0; ti<2; ti++){
        Eigen::Vector3i newTrig;
        for (int tv = 0; tv<3; tv++){
          int vi = em->e[ii]->at(HexFaces[fi][trigs[ti][tv]]);
          newTrig[tv] = e2surf[vi];
        }
        surf->t.push_back(newTrig);
        surf->ei.push_back(ii);
        surf->fi.push_back(fi);
      }
    }
  }
}

void savePly(TrigMesh * tm, std::string filename)
{
  FileUtilOut out(filename);
  out.out << "ply\nformat ascii 1.0\nelement vertex " << tm->v.size() << "\n";
  out.out << "property float x\nproperty float y\nproperty float z\n";
  out.out << "property uchar red\nproperty uchar green\nproperty uchar blue\n";
  out.out << "element face " << tm->t.size() << "\n";
  out.out << "property list uchar int vertex_index\nend_header\n";
  //compute average vertex color;
  std::vector<Eigen::Vector3f> vcolor(tm->v.size(), Eigen::Vector3f::Zero());
  std::vector<int> vcnt(tm->v.size(),0);
  for (size_t i = 0; i < tm->t.size(); i++){
    for (int j = 0; j < 3; j++){
      int vidx = tm->t[i][j];
      vcolor[vidx] += tm->tcolor[i].row(j);
      vcnt[vidx] ++;
    }
  }
  for (size_t i = 0; i < tm->v.size(); i++){
    vcolor[i] = (1.0 / vcnt[i]) * vcolor[i];
    out.out << tm->v[i][0] << " " << tm->v[i][1] << " " << tm->v[i][2]<<" ";
    out.out << (unsigned int)(vcolor[i][0] * 255) << " " << (unsigned int)(vcolor[i][1] * 255) << " " << (unsigned int)(vcolor[i][2] * 255) << "\n";
  }
  for (size_t i = 0; i < tm->t.size(); i++){
    out.out << "3 " << tm->t[i][0] << " " << tm->t[i][1] << " " << tm->t[i][2] << "\n";
  }
  out.close();
}

void saveVelObj(ElementMesh * em, std::string filename)
{
  FileUtilOut out(filename);
  TrigMesh arrow;
  FileUtilIn arrowIn("arrow.obj");
  arrow.load(arrowIn.in);
  arrowIn.close();
  out.close();
  TrigMesh all;
  float maxV = 0;
  for (int i = 0; i < em->v.size(); i++){
    maxV = std::max(maxV, (float)em->v[i].norm());
  }
  float eleSize = em->X[1][2] - em->X[0][2];
  //float vscale = eleSize / maxV;
  //std::cout << "vscale " << vscale << "\n";
  float vscale = 0.005;
  std::vector<Eigen::Vector3d> v0 = arrow.v;
  for (int i = 0; i < em->v.size(); i++){
    arrow.v = v0;
    Eigen::Vector3d x(1, 0, 0);
    float len = em->v[i].norm();
    if (len > 1e-6){
      x = em->v[i]/len;
    }
    Eigen::Vector3d y = Eigen::Vector3d(0, 1, 0);
    if (x[1]>0.9){
      y = Eigen::Vector3d(1, 0, 0);
    }
    Eigen::Vector3d z = x.cross(y).normalized();
    y = z.cross(x).normalized();


    for (size_t j = 0; j < arrow.v.size(); j++){
      arrow.v[j] = vscale * len*( arrow.v[j][0] * x + arrow.v[j][1] * y + arrow.v[j][2] * z);
      arrow.v[j] += em->x[i];
    }
    all.append(arrow);
  }

  all.save_obj(filename.c_str());
}

Eigen::Vector3f colormapBW(double val)
{
  Eigen::Vector3f color = ((float)val) * Eigen::Vector3f::Ones();
  return color;
}

/// @param val assumed to be between 0 and 1.
Eigen::Vector3f colormap(double val)
{
  Eigen::Vector3f color = Eigen::Vector3f::Zero();
  //red
  if (val >= 0.35 && val <= 0.6){
    color[0] = 4 * (val - 0.35);
  }
  else if (val >0.6 && val <= 0.85){
    color[0] = 1;
  }
  else if(val >0.85){
    color[0] = 1 - 3 * (val - 0.85);
  }

  //green
  if (val >= 0.1 && val <= 0.35){
    color[1] = 4 * (val - 0.1);
  }
  else if (val >0.35 && val <= 0.6){
    color[1] = 1;
  }
  else if (val > 0.6 && val <= 0.85){
    color[1] = 1 - 4 * (val - 0.6);
  }

  //blue
  if (val < 0.1){
    color[2] = 0.6 + 4 * val;
  }
  else if (val >=0.1 && val <= 0.35){
    color[2] = 1;
  }
  else if (val > 0.35 && val <= 0.6){
    color[2] = 1 - 4 * (val - 0.35);
  }

  return color;
}

//linearly map values to [0 1].
void normalize(std::vector<std::vector<double> > & a)
{
  double mn= a[0][0], mx = a[0][0];
  for (int i = 0; i < (int)a.size(); i++){
    for (int j = 0; j < (int)a[i].size(); j++){
      if (a[i][j] > mx){
        mx = a[i][j];
      }
      if (mn>a[i][j]){
        mn = a[i][j];
      }
    }
  }
  double eps = 1e-20;
  mx += eps;
  ///@TODO experimental.
  //mx = 5.1e5;
  //mx = 0.0065;
  //mx = 0.15;
  //mx = 60000;
  //mx = 320000;
  //mx = 1e7;
  std::cout << "normalize mn mx " << mn << " " << mx << "\n";
  for (int i = 0; i < (int)a.size(); i++){
    for (int j = 0; j < (int)a[i].size(); j++){
      //a[i][j] = (a[i][j] - mn) / (mx-mn) ;
      a[i][j] = a[i][j]  / mx ;
    }
  }
}

void hexToTrigStress(ElementMesh * em, TrigMesh * tm,
  std::vector<std::vector<double> > & P)
{
  quadratureStress(em, P);
  normalize(P);
  tm->tcolor.resize(tm->t.size(), 0.5f*Eigen::Matrix3f::Ones());

  for (int i = 0; i < (int)tm->t.size(); i++){
    int eidx = tm->ei[i];
    std::vector<int> verts = em->e[eidx]->verts();
    for (int j = 0; j < 3; j++){
      int mvidx = tm->t[i][j];      
      //vertex index in element mesh.
      int evidx = tm->vidx[mvidx];
      int localIdx = find(evidx, verts);
      double val = P[eidx][localIdx];
      Eigen::Vector3f color = colormap(val);
      //std::cout << color << "\n";
      tm->tcolor[i].row(j) = color;
    }    
  }
}

std::vector<std::vector<Eigen::Vector3d> >
loadModes(const char * filename, int nvert, int dim)
{
  std::ifstream in(filename);
  if (!in.good()){
    std::cout << "Cannot open " << filename << "\n";
  }
  std::string line;
  getline(in, line);
  std::stringstream ss(line);
  std::vector<double> arr;
  double num;
  while (ss >> num){
    arr.push_back(num);
  }
  int nModes = arr.size();
  std::cout << nModes << " modes\n";
  std::vector<std::vector<Eigen::Vector3d> > modes(nModes, std::vector<Eigen::Vector3d>(nvert));
  for (int ii = 0; ii<nModes; ii++){
    modes[ii][0][0] = arr[ii];
  }
  for (int ii = 0; ii<nvert; ii++){
    for (int jj = 0; jj<dim; jj++){
      if (ii == 0 && jj == 0){
        continue;
      }
      for (int kk = 0; kk<nModes; kk++){
        in >> num;
        modes[kk][ii][jj] = num;
      }
    }
  }
  in.close();
  return modes;
}

void savePartObj(ElementMesh * em, int part, std::string filename)
{
  std::vector<Element*> e0 = em->e;
  std::vector<Element*> subset;
  for (size_t i = 0; i < em->e.size(); i++){
    if (em->me[i] == part){
      subset.push_back(em->e[i]);
    }
  }
  if (subset.size() == 0){
    em->e = e0;
    return;
  }
  em->e = subset;
  TrigMesh tm;
  hexToTrigMesh(em, &tm);
  tm.save_obj(filename.c_str());
  em->e = e0;
}

void saveSimState(const ElementMesh * em, std::ostream & out)
{
  save_arr(em->x, out);
  save_arr(em->v, out);
  save_arr(em->fe, out);
  save_arr(em->fc, out);
}

void loadSimState(ElementMesh * em, std::istream & in, bool loadDisplacement)
{
  load_arr(em->x, in);
  if (loadDisplacement){
    for (size_t i = 0; i < em->x.size(); i++){
      em->x[i] += em->X[i];
    }
  }
  load_arr(em->v, in);
  load_arr(em->fe, in);
  em->fc.resize(em->fe.size());
  load_arr(em->fc, in);
}

double vonMises(Eigen::Matrix3d sigma)
{
  double pi = (1.0/3.0) * sigma.trace();
  Eigen::Matrix3d s = sigma - pi * Eigen::Matrix3d::Identity();
  double J2 = 0.5 * (s*s.transpose()).trace();
  return std::sqrt(J2);
}

void quadratureStress(ElementMesh * em, std::vector<std::vector<double> > & P)
{
  P.resize(em->e.size());
  int nQ = (int)em->e[0]->detJ.size();
  for (size_t i = 0; i < em->e.size(); i++){
    int matIdx = em->me[i];
    Material * m = em->m[matIdx];
    P[i].resize(nQ);
    for (int j = 0; j < nQ; j++){
      Eigen::Matrix3d pk1 = m->getStress(j, i, em);
      Eigen::Matrix3d F = em->e[i]->defGradCached(j, em->X, em->x);
      double J = F.determinant();
      if (J <= 0){
        //inverted
        P[i][j] = 0;
        continue;
      }
      Eigen::Matrix3d sigma = (1.0/J)*F*pk1;
      MaterialQuad *mq = (MaterialQuad *)m;
      
      double v = 0;
      for (int k = 0; k < 3; k++){
        v += std::abs((F*F.transpose()).trace() - 3);
      }
      //v += std::max(0.0, F.determinant()-1);
      v = vonMises(sigma);
      P[i][j] = v;
    }
  }
}

void saveAbq(ElementMesh * em, std::string filename, float scale)
{
  FileUtilOut out(filename);
  out.out << "*NODE\n";
  for (unsigned int ii = 0; ii<em->X.size(); ii++){
    Eigen::Vector3d v = em->X[ii];
    v *= scale;
    out.out << (ii + 1) << ", " << v[0] << ", " << v[1] << ", " << v[2] << "\n";
  }
  out.out << "*ELEMENT, TYPE=C3D8\n";
  int abqOrder[8] = { 0, 2, 6, 4, 1, 3, 7, 5 };
  for (unsigned int ii = 0; ii<em->e.size(); ii++){
    out.out << (ii + 1);
    for (int jj = 0; jj<em->e[ii]->nV(); jj++){
      out.out << ", " << (em->e[ii]->at(abqOrder[jj]) + 1);
    }
    out.out << "\n";
  }
  out.out << "*ELSET,ELSET=EALL,GENERATE\n";
  out.out << "  1," << (em->e.size()) << "\n";
  out.close();
}

///\brief normalize so that largest displacement is 10mm.
void normalize(std::vector<Eigen::Vector3d> & v, int keyIdx)
{
  double mx = v[0].squaredNorm();
  int midx = 0;
  for (size_t i = 1; i < v.size(); i++){
    double n = v[i].squaredNorm();
    if (n > mx){
      mx = n;
      midx = i;
    }
  }
  double norm = std::sqrt(mx);
  //double norm = v[keyIdx].norm();
  double sum = 0;
  std::cout << v[midx] << "\n";
  for (size_t i = 0; i < v.size(); i++){
    
    for (int j = 0; j < 3; j++){
      sum += v[i][j];
    }
  }
  mx = 0;
  
  if (sum < 0){
    norm = -1 * norm;
  }
  std::cout << "sum " << sum << "\n";
  for (size_t i = 0; i < v.size(); i++){
    v[i] /= 130*norm;
  }
}

int topVert(const std::vector < Eigen::Vector3d> & X)
{
  int idx;
  double max = 0;
  for (size_t i = 0; i < X.size(); i++){
    if (X[i][1] > max){
      max = X[i][1];
      idx = i;
    }
  }
  return idx;
}

void saveModalMeshes(ElementMesh * em, std::string modefile)
{
  std::vector<std::vector<Eigen::Vector3d> > modes = loadModes(modefile.c_str(), em->x.size(), em->dim);

  std::cout << "modes " << modes.size() << "\n";
  //FileUtilIn in("0/0m33.txt");
  //load_arr(em->x, in.in);
  //in.close();

  std::vector<Eigen::Vector3d> u0(em->X.size());
  for (size_t i = 0; i < em->x.size(); i++){
    u0[i] = em->x[i] - em->X[i];
  }
  em->x = em->X;
  TrigMesh tm;
  hexToTrigMesh(em, &tm);
  int keyIdx = 0;
  keyIdx = topVert(em->X);
  std::cout << "keyIdx" << keyIdx << "\n"<< em->X[keyIdx] << "\n";
  for (int modeNum = 0; modeNum < modes.size(); modeNum++){
    em->x = em->X;
    std::cout << "mode " << modeNum << "\n";
    normalize(modes[modeNum], keyIdx);
    double mag = 1;
    if (modeNum == 7){
      mag = -1;
    }
    // double mag = dot(modes[modeNum], u0);
    //mag = 0.6 * mag / std::abs(mag);
    for (size_t i = 0; i < em->x.size(); i++){
      //em->x[i] += mag * modes[modeNum][i];
      em->x[i] += 3* modes[modeNum][i];
    }

    updateSurfVert(em->x, &tm, tm.vidx);
    std::string filename = "mode" + std::to_string(modeNum) + "-0.obj";
    savePartObj(em, 0, filename);
    //filename = "mode" + std::to_string(modeNum) + "-1.obj";
    //savePartObj(em, 1, filename);
  }

  ////TrigMesh tm;
  ////hexToTrigMesh(em, &tm);
  ////updateSurfVert(em->x, &tm, tm.vidx);
  ////hexToTrigStress(em, &tm);

  ////savePly(&tm, "volume.ply");
  ////savePartObj(em, 0, "modeProj-0.obj");
  ////savePartObj(em, 1, "modeProj-1.obj");
}

void saveStiffness(ElementMesh * em)
{
  Eigen::SparseMatrix <double> K = em->getStiffnessSparse();
  Eigen::SparseMatrix <double> M = em->M;
  zeroOffDiag(K, em->fixedDof);
  zeroOffDiag(M, em->fixedDof);
  std::ofstream out("K.txt");
  write_vega_lists(out, K);
  out.close();
  out.open("M.txt");
  write_vega_lists(out, M);
  out.close();
}

///@brief moves rest pose X to origin
///and sets x=X.
void placeAtOrigin(ElementMesh * em)
{
  Eigen::Vector3d mn, mx;
  BBox(em->X, mn, mx);
  mn[0] = mx[0];
  for (size_t i = 0; i < em->x.size(); i++){
    em->X[i] -= mn;
  }
  em->x = em->X;
}

void scaleMesh(ElementMesh * em, float scale)
{
  for (unsigned int ii = 0; ii<em->X.size(); ii++){
    em->X[ii] *= scale;
  }
  em->x = em->X;
  for (unsigned int ii = 0; ii<em->depth.size(); ii++){
    em->depth[ii] *= scale;
  }
}

void scaleMesh(ElementMesh * em, const std::vector<float> & scale)
{
  for (unsigned int ii = 0; ii<em->X.size(); ii++){
    for (int j = 0; j < em->dim; j++){
      em->X[ii][j] *= scale[j];
    }
  }
  em->x = em->X;
}

void scaleMeshx(ElementMesh * em, const std::vector<float> & scale)
{
  for (unsigned int ii = 0; ii<em->X.size(); ii++){
    for (int j = 0; j < em->dim; j++){
      em->x[ii][j] *= scale[j];
    }
  }
}

void translate(ElementMesh * em, const Eigen::Vector3d & vec)
{
  for (size_t i = 0;i< em->x.size(); i++){
    em->x[i] += vec;
  }
}

///@param O new corner (0,0,0) of bounding box.
///@param Frame. Each column is an axis in world space.
void applyFrame(ElementMesh * em, const Eigen::Vector3d & O,
  const Eigen::Matrix3d & Frame)
{
  Eigen::Vector3d mn, mx;
  BBox(em->x, mn, mx);
  for (size_t i = 0; i < em->x.size(); i++){
    em->x[i] += O - mn;
    em->x[i] = Frame * (em->x[i] - O) + O;
  }
}

bool hitWall(float wallDist, const std::vector<Eigen::Vector3d> & x, int dim, int sign)
{
  bool intersect = false;
  for (size_t i = 0; i < x.size(); i++){
    if (sign*(wallDist - x[i][dim]) < 0){
      intersect = true;
      break;
    }
  }
  return intersect;
}
