#include "ElementMesh.hpp"
#include "ElementHex.hpp"
#include "ElementQuad.hpp"
#include "ElementMeshUtil.hpp"
#include "Material.hpp"
#include "Element.hpp"
#include "Quadrature.hpp"
#include "ArrayUtil.hpp"
#include "TrigMesh.hpp"

#include "femError.hpp"
#include <map>
#include <iostream>
#include <fstream>

#define PI 3.14159265

typedef Eigen::Triplet<double> Tripletd;

void
ElementMesh::grid2D(std::vector<std::vector< int> > & grid)const
{
  Eigen::Vector3d mn, mx;
  BBox(X, mn, mx);
//  std::cout<<"Bounding box "<<mn[0] << " " << mn[1] << " " <<mn[2]<<"\n";
//  std::cout<< mx[0] << " " << mx[1] << " " <<mx[2]<<"\n";
  double gridlen = X[ e[0]->at(1) ][1] - X[ e[0]->at(0) ][1];
//  std::cout<<"Grid len "<<gridlen<<"\n";

  int nx, ny;
  nx = (int)( (mx[0] - mn[0])/gridlen + 0.5 );
  ny = (int)( (mx[1] - mn[1])/gridlen + 0.5 );

  grid.resize(nx);
  for(unsigned int ii = 0; ii<grid.size(); ii++){
    grid[ii].resize(ny, -1);
  }

  for(unsigned int ii = 0; ii<e.size(); ii++){
    int vidx = e[ii]->at(0);
    Eigen::Vector3d v = X[vidx];
    int xx = (int)( (v[0]-mn[0]) / gridlen + 0.5 );
    int yy = (int)( (v[1]-mn[1]) / gridlen + 0.5 );
    grid[xx][yy] = ii;
  }
}

void
ElementMesh::grid3D(std::vector<std::vector< std::vector<int> > > & grid)const
{
  Eigen::Vector3d mn, mx;
  BBox(X, mn, mx);
//  std::cout<<"Bounding box "<<mn[0] << " " << mn[1] << " " <<mn[2]<<"\n";
//  std::cout<< mx[0] << " " << mx[1] << " " <<mx[2]<<"\n";
  double gridlen = X[ e[0]->at(2) ][1] - X[ e[0]->at(0) ][1];
//  std::cout<<"Grid len "<<gridlen<<"\n";

  int nx, ny, nz;
  nx = (int)( (mx[0] - mn[0])/gridlen + 0.5 );
  ny = (int)( (mx[1] - mn[1])/gridlen + 0.5 );
  nz = (int)( (mx[2] - mn[2])/gridlen + 0.5 );

  grid.resize(nx);
  for(int ii = 0; ii<nx; ii++){
    grid[ii].resize(ny);
    for(int jj = 0; jj<ny; jj++){
      grid[ii][jj].resize(nz, -1);
    }
  }

  for(unsigned int ii = 0; ii<e.size(); ii++){
    int vidx = e[ii]->at(0);
    Eigen::Vector3d v = X[vidx];
    int xx = (int)( (v[0]-mn[0]) / gridlen + 0.5 );
    int yy = (int)( (v[1]-mn[1]) / gridlen + 0.5 );
    int zz = (int)( (v[2]-mn[2]) / gridlen + 0.5 );
    grid[xx][yy][zz] = ii;
  }
}

///@TODO only loads hex mesh
void
ElementMesh::load(std::istream & in, float scale)
{
  int nEle , nVert;
  std::string token;
  std::cout<<"load mesh...\n";
  in>>token;
  in>>nVert;
  in>>token;
  in>>nEle;

  X.resize(nVert);
  e.resize(nEle);
  std::cout<<nVert<<" vertices\n"<<nEle<<" elements.\n";
  for(int ii =0 ; ii<nVert; ii++){
    in>>X[ii][0]>>X[ii][1]>>X[ii][2];
  }

  for(int ii =0; ii<nEle; ii++){
    int nV;
    in>>nV;
    Element * ele = 0;
    if(nV==8){
      ele = new ElementHex();
    }else if(nV==4){
      ele = new ElementQuad();
    }
    e[ii] = ele;
    for(int jj =0 ; jj<nV; jj++){
      in>>(*ele)[jj];
    }
  }

  initArrays();
  place(Eigen::Vector3d::Zero());
  //read optional material specification;
  in>>token;
  in>>nEle;
  int nPart = 1;
  if(token == "parts"){
    for(int ii =0; ii<nEle; ii++){
      int p;
      in>>p;
      if(p>nPart){
        nPart = p;
      }
      me[ii] = p - 1;
    }

    in>>token;
    in>>nEle;
    for(int ii =0; ii<nEle; ii++){
      float d;
      in>>d;
      depth[ii] = d;
    }
  }

  materialColor.resize(nPart);
  materialColor[0] = Eigen::Vector3f(0.7, 0.7, 0.9);
  if (nPart > 1){
    materialColor[1] = Eigen::Vector3f(0, 0, 0);
  }
  scaleMesh(this, scale);
  //in 3d meshes depth is set to 1.
  if(dim==3){
    std::fill(depth.begin(), depth.end(), 1.0);
  }
}

void ElementMesh::place(Eigen::Vector3d origin)
{
  Eigen::Vector3d minPos(10, 10, 10);
  if(dim==2){
    minPos[2] = 0;
  }
  for(unsigned ii =0 ; ii<X.size(); ii++){
    for(int jj =0 ; jj<dim; jj++){
      if(minPos[jj]>X[ii][jj]){
        minPos[jj] = X[ii][jj];
      }
    }
  }
  for(auto it = X.begin(); it!=X.end(); it++){
    (*it) -= (minPos-origin);
  }
  x = X;
}

void ElementMesh::saveMesh(std::ostream &out)
{
  out<<"#verts "<<X.size();
  out<<"\n#elts "<<e.size()<<"\n";
  for(unsigned int ii =0; ii<X.size(); ii++){
    out<<X[ii][0]<<" "<<X[ii][1]<<" "<<X[ii][2]<<"\n";
  }
  for(unsigned int ii =0; ii<e.size(); ii++){
    out<<e[ii]->nV();
    for(int jj =0 ; jj<e[ii]->nV(); jj++){
      out<<" "<<e[ii]->at(jj);
    }
    out<<"\n";
  }

  out<<"parts\n"<<e.size()<<"\n";
  for(unsigned int ii = 0; ii<e.size(); ii++){
    out<<me[ii]+1<<"\n";
  }
  out<<"depth\n"<<e.size()<<"\n";
  for(unsigned int ii = 0; ii<e.size(); ii++){
    out<<depth[ii]<<"\n";
  }

}

void ElementMesh::initArrays()
{
  x  = X;
  v  = std::vector<Eigen::Vector3d>(X.size(), Eigen::Vector3d (0,0,0));
  fe = std::vector<Eigen::Vector3d>(X.size(), Eigen::Vector3d (0,0,0));
  lb = std::vector<Eigen::Vector3d>(X.size(), Eigen::Vector3d (-10, 0,-10));
  fixedDof = std::vector<int>(X.size()*dim);
  me = std::vector<int>(e.size(), 0);
  depth = std::vector<double>(e.size(), 1);
  materialColor = std::vector<Eigen::Vector3f>(1, Eigen::Vector3f(0,0,1));
}

void ElementMesh::initElements(const Quadrature *q)
{
  Element * ele = e[0];
  //pre-compute shape function gradients
  if (Element::gradN.size() != q->x.size()){
    Element::gradN.resize(q->x.size());
    for (int qq = 0; qq < q->x.size(); qq++){
      Element::gradN[qq].resize(ele->nV());
      for (int ii = 0; ii < ele->nV(); ii++){
        Element::gradN[qq][ii] = ele->shapeFunGrad(ii, q->x[qq]);
      }
    }
  }

  //compute per-element jacobian that depends on the reset shape X
  for(unsigned int ii = 0; ii<e.size(); ii++){
    ele = e[ii];
    ele->Jinv.resize(q->x.size());
    ele->detJ.resize(q->x.size());
    ele->JdN.resize(q->x.size());
    for (int qq = 0; qq < q->x.size(); qq++){
      Eigen::Matrix3d J = ele->getJ(qq, X);
      if(dim==2){
        J(2,2) = 1;
      }
      ele->Jinv[qq] = J.inverse();
      ele->detJ[qq] = J.determinant();
      ele->JdN[qq].resize(ele->nV());
      for(int jj = 0; jj<ele->nV(); jj++){
        ele->JdN[qq][jj] = ele->Jinv[qq].transpose() * ele->gradN[qq][jj];
      }
    }
  }

  computeMass();
  fc.resize(X.size(),Eigen::Vector3d(0,0,0));
}

void ElementMesh::addMaterial(Material*_m)
{
  m.push_back(_m);
  _m->init(this);
}

int ElementMesh::check()
{
  if(fe.size()!=x.size()){
    std::cout<<"external force and dof size differ\n";
    return -1;
  }
  if(e.size()!=me.size()){
    std::cout<<"material assignment and num elements differ\n";
    return -1;
  }
  for(unsigned int ii = 0;ii<me.size();ii++){
    if(me[ii]<0 || me[ii]>=m.size()){
      std::cout<<"Material index out of bounds\n";
      std::cout<<"ele: "<<ii<<", mat: "<<me[ii]<<"\n";
      return -1;
    }
  }
    
  if(fixedDof.size() != dim * x.size()){
    std::cout<<"fixed array differ in size to vertex array \n";
    return -1;
  }
  return 0;
}

double ElementMesh::getAllEnergy(bool addgravity)
{
  double ene =getEnergy(addgravity);
  ene += getKineticEnergy();
  return ene;
}

double ElementMesh::getEnergyElastic()
{
  double ene = 0;
  for (unsigned int ii = 0; ii<e.size(); ii++){
    ene += getEnergyEle(ii);
  }
  return ene;
}

double ElementMesh::getEnergy(bool addgravity)
{
  double ene = 0;
  double Ei, Ee;
  for(unsigned int ii = 0;ii<e.size();ii++){
    ene += getEnergyEle(ii);
    if(fem_error){
      return -1;
    }
  }
  Ei = ene;
  
  //energy from external forces
  for(unsigned int ii = 0;ii<fe.size();ii++){
    ene -= fe[ii].transpose()*(x[ii] - X[ii]);
  }
  Ee = Ei - ene;
  //std::cout << "E " << Ei <<" " << Ee << "\n";
  if(addgravity){
    for(unsigned int ii = 0;ii<fe.size();ii++){
      ene -= mass[ii] * x[ii].transpose() * G;
    }
  }
  return ene;
}

double ElementMesh::getKineticEnergy(bool addgravity)
{
  Eigen::VectorXd vec(dim * v.size());
  for(unsigned int ii = 0; ii<v.size(); ii++){
    for(int jj = 0; jj<dim; jj++){
      vec[ii*dim + jj] = v[ii][jj];
    }
  }
  double ene = 0.5*vec.transpose()*M*vec;
  if (addgravity){
    for (unsigned int ii = 0; ii<fe.size(); ii++){
      ene -= mass[ii] * x[ii].transpose() * G;
    }
  }
  return ene;
}

double ElementMesh::getEnergyEle(int eIdx)
{
  return m[me[eIdx]]->getEnergy(eIdx, this);
}

std::vector<Eigen::Vector3d> ElementMesh::getForceEle(int eIdx)
{
  return m[me[eIdx]]->getForce(eIdx, this);
}

std::vector<Eigen::Vector3d> ElementMesh::getForce(bool add_gravity)
{
  std::vector<Eigen::Vector3d> force(x.size(), Eigen::Vector3d::Zero());
  for(unsigned int ii = 0;ii<e.size();ii++){
    std::vector<Eigen::Vector3d> fele = getForceEle(ii);
    if(fem_error){
      return force;
    }
    for(int jj = 0; jj<e[ii]->nV(); jj++){
      force[ e[ii]->at(jj) ] += fele[jj];
    }
  }
  for(unsigned int ii= 0;ii<fe.size();ii++){
    force[ii] += fe[ii];
  }

  if(add_gravity){
    for(unsigned int ii= 0;ii<fe.size();ii++){
      force[ii] += mass[ii] * G;
    }
  }
  return force;
}

Eigen::VectorXd
ElementMesh::getdForceEle(int eIdx, const Eigen::VectorXd & dx)
{
///@TODO unimplemented
  Eigen::VectorXd dfe;
  return dfe;
}

Eigen::VectorXd
ElementMesh::getdForce(const Eigen::VectorXd & dx)
{
  Eigen::VectorXd dfe;
  Eigen::VectorXd df=Eigen::VectorXd::Zero(dx.size());
  for(unsigned int ei = 0; ei<e.size(); ei++){
    dfe = getdForceEle(ei, dx);
    Element * ele = e[ei];
    for(int eleVi = 0; eleVi<ele->nV(); eleVi++){
      int globalVi = ele->at(eleVi);
      df.block(dim * globalVi,0, dim,1) += dfe.block(dim*eleVi,0,dim,1);
    }
  }
  return df;
}



///@TODO: only works for rectangular elements
void ElementMesh::computeMass()
{
//  std::cout<<"Compute mass...\n";

  totalMass =0;
  totalVol = 0;

  int N = dim * (int)X.size();
  M = Eigen::SparseMatrix<double>(N,N);
  M.reserve(Eigen::VectorXi::Constant(N, 81));
  mass= std::vector<float>(X.size(),0);
  std::vector<std::map<int, double> > massMap(N);

  for(unsigned int ii = 0;ii<e.size();ii++){
    Element * ele = e[ii];
    double vol = e[ii]->getVol(X);
    int mat = me[ii];
    if (mat >= m.size()){
      mat = m.size() - 1;
    }
    double eleMass = vol*m[mat]->density * depth[ii];
    totalVol  += (float)(vol*depth[ii]);
    totalMass += (float)(eleMass);
    Eigen::MatrixXd M0 = m[mat]->getMassMatrix(ii, this);

    for(int jj = 0; jj<ele->nV(); jj++){
      int vj = ele->at(jj);
      mass[vj] += (float)(eleMass/ele->nV());
      for(int kk = 0; kk<ele->nV(); kk++){
        int vk = ele->at(kk);
        for(int dim1= 0 ;dim1<dim; dim1++){
          for(int dim2= 0 ;dim2<dim; dim2++){
            double val = M0(dim * jj + dim1, dim*kk + dim2);
            val *= m[mat]->density * depth[ii];
            M.coeffRef(dim*vj + dim1, dim*vk + dim2) += val;
          }
        }
      }
    }
  }
  M.makeCompressed();
}

Eigen::Vector3d
ElementMesh::centerOfMass()const
{
  Eigen::Vector3d center(0,0,0);
  double totalMass = 0;
  for(unsigned int ii = 0; ii<x.size(); ii++){
    center += x[ii] * mass[ii];
    totalMass += mass[ii];
  }
  center = (1.0/totalMass) * center;
  return center;
}

Eigen::Vector3d
ElementMesh::angularVelocity(Eigen::Vector3d center)const
{
  Eigen::Vector3d omega(0,0,0);
  double mtotal = 0;
  for(unsigned int ii = 0; ii<x.size(); ii++){
    Eigen::Vector3d r = x[ii] - center;
    double len = r.squaredNorm();
    omega += mass[ii]*(x[ii]-center).cross(v[ii])/len;
    //weight by vertex mass
    mtotal += mass[ii];
  }
  return omega/mtotal;
}

Eigen::Matrix3d
ElementMesh::inertiaTensor(Eigen::Vector3d center)const
{
  Eigen::Matrix3d Inertia=Eigen::Matrix3d::Zero();
  Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
  const Quadrature * quad = &Quadrature::Gauss2;
  for(unsigned int ei = 0; ei<e.size(); ei++){
    Element * ele = e[ei];
    //3d quadrature rules sum to 8.
    double rho = m[me[ei]]->density * ele->getVol(X) / 8;
    for(unsigned int qi = 0; qi<quad->x.size(); qi++){
      std::vector<double> N = ele->shapeFun(quad->x[qi]);
      Eigen::Vector3d xq(0,0,0);
      for(unsigned int ii = 0; ii<N.size(); ii++){
        xq += N[ii] * x[ele->at(ii)];
      }
      xq -= center;
      Inertia += quad->w[qi] * rho * (xq.dot(xq) * I - xq * xq.transpose());
    }
  }
  return Inertia;
}

Eigen::Vector3d
ElementMesh::angularMomentum(Eigen::Vector3d center)const
{
  Eigen::Vector3d L(0,0,0);

  Eigen::VectorXd vec(dim * v.size());
  for(unsigned int ii=0; ii<v.size(); ii++){
    for(int jj = 0; jj<dim; jj++){
      vec[ii*dim + jj] = v[ii][jj];
    }
  }
  Eigen::VectorXd pvec = M*vec;
  Eigen::Vector3d p=Eigen::Vector3d::Zero();
  for(unsigned int ii = 0; ii<x.size(); ii++){
    for(int jj =0; jj<dim; jj++){
      p[jj] = pvec[ii*dim + jj];
    }
    L += (x[ii]-center).cross(p);
  }

  return L;
}

double
ElementMesh::getTopAng()
{
  if(topVerts.size()<2){
      return 0;
  }
//  std::cout<<x[topVerts[0]][0]<<" "<<x[topVerts[0]][1]<<" "<<x[topVerts[0]][2]<<"\n";
//  std::cout<<x[topVerts[1]][0]<<" "<<x[topVerts[1]][1]<<" "<<x[topVerts[1]][2]<<"\n";
  Eigen::Vector3d vec = x[topVerts[0]] - x[topVerts[1]];
  double ang = std::atan2(vec[1], vec[0]);
  return ang;
}

Eigen::Vector3d
ElementMesh::linearVelocity()const
{
  Eigen::Vector3d p(0,0,0);
  double mtotal = 0;
  for(unsigned int ii =0 ; ii<v.size(); ii++){
    p += mass[ii] * v[ii];
    mtotal += mass[ii];
  }
  return p/mtotal;
}


Eigen::Vector3d
ElementMesh::linearMomentum()const
{
  Eigen::VectorXd vec(dim * v.size());
  for(unsigned int ii=0; ii<v.size(); ii++){
    for(int jj = 0; jj<dim; jj++){
      vec[ii*dim + jj] = v[ii][jj];
    }
  }
  Eigen::VectorXd pvec = M*vec;
  Eigen::Vector3d p=Eigen::Vector3d::Zero();

  for(unsigned int ii =0 ; ii<v.size(); ii++){
    for(int jj = 0; jj<dim; jj++){
      p[jj] += pvec[ii*dim+jj];
    }
  }
  return p;
}

ElementMesh::ElementMesh():dim(3),
  G(Eigen::Vector3d(0,-9.8,0)),
  dampAlpha(1e-1f), dampBeta(0),
  specialVert(-1),
  forceDrawingScale(1e0f),
  depthColorScale(30),u(0)
{}

ElementMesh::~ElementMesh()
{
  for(unsigned int ii = 0; ii<e.size();ii++){
    delete e[ii];
  }
}
