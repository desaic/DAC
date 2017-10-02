#include "TrigMesh.hpp"
#include <fstream>
#include <iostream>
#include <algorithm>
#ifdef _WIN32
#define NOMINMAX //Stop errors with std::max
#include <windows.h>
#endif
#include <stdio.h>
#include <cstdlib>
#include <utility>
#include <map>
#include <sstream>
#include <string.h>
#include "util.h"
//#include "ElementMesh.hpp"

void makeCube(TrigMesh & m, const Vector3s & mn,
    const Vector3s mx)
{
  Vector3s ss = mx -mn;
  Matrix33sr transform = Matrix33sr::Identity();
  for(int ii =0 ;ii<3;ii++){
    transform(ii,ii) = ss[ii];
  }
  m=UNIT_CUBE;
  for(unsigned int ii = 0;ii<m.v.size();ii++){
    m.v[ii] = mn + transform * m.v[ii];
  }
}

void TrigMesh::append(const TrigMesh & m)
{
  unsigned int offset = v.size();
  unsigned int ot = t.size();
  v.insert(v.end(),m.v.begin(),m.v.end());
  t.insert(t.end(),m.t.begin(), m.t.end());
  for(unsigned int ii = ot;ii<t.size();ii++){
    for(int jj = 0 ;jj<3;jj++){
      t[ii][jj] += offset;
    }
  }
}

TrigMesh & TrigMesh::operator= (const TrigMesh& m)
{
  v = m.v;
  t = m.t;
  n = m.n;
  return *this;
}

///@brief cube [0,1]^3
Vector3s CUBE_VERT[8]={
    Vector3s (0, 0, 0),
    Vector3s (1, 0, 0),
    Vector3s (1, 1, 0),
    Vector3s (0, 1, 0),
    Vector3s (0, 0, 1),
    Vector3s (1, 0, 1),
    Vector3s (1, 1, 1),
    Vector3s (0, 1, 1)
};

Vector3i CUBE_TRIG[12]={Vector3i(0,3,1),
    Vector3i(1, 3, 2),
    Vector3i(5, 4, 0),
    Vector3i(5, 0, 1),
    Vector3i(6, 5, 1),
    Vector3i(1, 2, 6),
    Vector3i(3, 6, 2),
    Vector3i(3, 7, 6),
    Vector3i(4, 3, 0),
    Vector3i(4, 7, 3),
    Vector3i(7, 4, 5),
    Vector3i(7, 5, 6)};

TrigMesh UNIT_CUBE(CUBE_VERT,CUBE_TRIG);

TrigMesh::TrigMesh():v(0),t(0){}

TrigMesh::TrigMesh(const std::vector<Vector3s> &_v,
    const std::vector<Vector3i> &_t):v(_v),t(_t)
{
  compute_norm();
}

TrigMesh::TrigMesh(const Vector3s * _v,
  const Vector3i * _t)
{
  v.assign(_v,_v+8);
  t.assign(_t,_t+12);
  
  compute_norm();
}

void TrigMesh::save(std::ostream & out, std::vector<Vector3s> * vert)
{
  std::string vTok("v");
  std::string fTok("f");
  std::string texTok("vt");
  char bslash='/';
  std::string tok;
  if(vert==0){
    vert = &v;
  }
  for(size_t ii=0;ii<vert->size();ii++){
    out<<vTok<<" "<<(*vert)[ii][0]<<" "<<(*vert)[ii][1]<<" "<<(*vert)[ii][2]<<"\n";
  }
  if(tex.size()>0){
    for(size_t ii=0;ii<tex.size();ii++){
      out<<texTok<<" "<<tex[ii][0]<<" "<<tex[ii][1]<<"\n";
    }
    for(size_t ii=0;ii<t.size();ii++){
      out<<fTok<<" "<<t[ii][0]+1<<bslash<<texId[ii][0]+1<<" "
      <<t[ii][1]+1<<bslash<<texId[ii][1]+1<<" "
      <<t[ii][2]+1<<bslash<<texId[ii][2]+1<<"\n";
    }
  }else{
    for(size_t ii=0;ii<t.size();ii++){
      out<<fTok<<" "<<t[ii][0]+1<<" "<<
          t[ii][1]+1<<" "<<t[ii][2]+1<<"\n";
    }
  }
  out<<"#end\n";
}

void TrigMesh::save(const char * filename)
{
  std::ofstream out;
  out.open(filename);
  save(out);
  out.close();
}


void TrigMesh::load(std::istream &in)
{
  read_obj(in);
}

void TrigMesh::read_obj(std::istream & f)
{
  std::string line;
  std::string vTok("v");
  std::string fTok("f");
  std::string texTok("vt");
  char bslash='/',space=' ';
  std::string tok;
  while(1) {
    std::getline(f,line);
    if(f.eof()) {
      break;
    }
    if(line == "#end"){
      break;
    }
    if(line.size()<3) {
      continue;
    }
    if(line.at(0)=='#') {
      continue;
    }
    std::stringstream ss(line);
    ss>>tok;
    if(tok==vTok) {
      Vector3s vec;
      ss>>vec[0]>>vec[1]>>vec[2];
      v.push_back(vec);
    } else if(tok==fTok) {
      bool hasTexture = false;
      if (line.find(bslash) != std::string::npos) {
        std::replace(line.begin(), line.end(), bslash, space);
        hasTexture = true;
      }
      std::stringstream facess(line);
      facess>>tok;
      std::vector<int> vidx;
      std::vector<int> texIdx;
      int x;
      while(facess>>x){
        vidx.push_back(x);
        if(hasTexture){
          facess>>x;
          texIdx.push_back(x);
        }
      }
      texIdx.resize(vidx.size());
      for(int ii = 0;ii<vidx.size()-2;ii++){
        Vector3i trig, textureId;
        trig[0] = vidx[0]-1;
        textureId[0] = texIdx[0]-1;
        for (int jj = 1; jj < 3; jj++) {
          trig[jj] = vidx[ii+jj]-1;
          textureId[jj] = texIdx[ii+jj]-1;
        }
        t.push_back(trig);
        texId.push_back(textureId);
      }
    } else if(tok==texTok) {
      Vector2s texcoord;
      ss>>texcoord[0];
      ss>>texcoord[1];
      tex.push_back(texcoord);
    }
  }
  std::cout<<"Num Triangles: "<< t.size()<<"\n";
}

void TrigMesh::read_ply(std::istream & f)
{
  std::string line;
  std::string vertLine("element vertex");
  std::string faceLine("element face");
  std::string texLine("property float s");
  std::string endHeaderLine("end_header");
  while(true) {
    std::getline(f,line);
    if(std::string::npos!=line.find(vertLine)) {
      break;
    }
  }
  std::string token;
  std::stringstream ss(line);
  ss>>token>>token;
  int nvert;
  ss>>nvert;
  bool hasTex=false;
  while(true) {
    std::getline(f,line);
    if(std::string::npos!=line.find(faceLine)) {
      break;
    }
    if(std::string::npos!=line.find(texLine)) {
      hasTex=true;
    }
  }
  std::stringstream ss1(line);
  ss1>>token>>token;
  int nface;
  ss1>>nface;
  while(true) {
    std::getline(f,line);
    if(std::string::npos!=line.find(endHeaderLine)) {
      break;
    }
  }

  v.resize(nvert);
  t.resize(nface);
  if(hasTex) {
    tex.resize(nvert);
  }
  for (int ii =0; ii<nvert; ii++) {
    for (int jj=0; jj<3; jj++) {
      f>>v[ii][jj];
    }
    if(hasTex) {
      for (int jj=0; jj<2; jj++) {
        f>>tex[ii][jj];
      }
      tex[ii][1]=1-tex[ii][1];;
    }
  }
  for (int ii =0; ii<nface; ii++) {
    int nidx;
    f>>nidx;
    for (int jj=0; jj<3; jj++) {
      f>>t[ii][jj];
    }
  }
}

//void TrigMesh::save_obj(const char * filename, int mat, ElementMesh * em)
//{
//  std::ofstream out(filename);
//  if (!out.good()){
//    std::cout << "cannot open output file" << filename << "\n";
//    return;
//  }
//  std::string vTok("v");
//  std::string fTok("f");
//  std::string texTok("vt");
//  char bslash = '/';
//  std::string tok;
//
//  for (size_t ii = 0; ii<v.size(); ii++){
//    out << vTok << " " << v[ii][0] << " " << v[ii][1] << " " << v[ii][2] << "\n";
//  }
//
//  for (size_t ii = 0; ii < t.size(); ii++){
//    int eidx = ei[ii];
//    int midx = em->me[eidx];
//    if (midx == mat){
//      out << fTok << " " << t[ii][0] + 1 << " " <<
//        t[ii][1] + 1 << " " << t[ii][2] + 1 << "\n";
//    }
//  }
//  out << "#end\n";
//  out.close();
//}

void TrigMesh::save_obj(const char * filename)
{
  std::ofstream out(filename);
  if(!out.good()){
    std::cout<<"cannot open output file"<<filename<<"\n";
    return;
  }
  save(out);
  out.close();
}

void TrigMesh::update()
{}

TrigMesh::TrigMesh(const char * filename,bool normalize)
{
  load_mesh(filename,normalize);
}


void TrigMesh::load_mesh(const char * filename, bool normalize)
{
  std::ifstream f ;
  f.open(filename);
  if(!f.good()) {
    std::cout<<"Error: cannot open mesh "<<filename<<"\n";
    return;
  }
  switch(filename[strlen(filename)-1]) {
  case 'y':
    read_ply(f);
    break;
  case 'j':
    read_obj(f);
    break;
  default:
    break;
  }
  if(normalize){
    rescale();
  }
  compute_norm();

  f.close();
}

void TrigMesh::rescale()
{
  if(v.size()==0){
    std::cout<<"empty mesh\n";
    return;
  }
  Vector3s mn=v[0],mx=v[0];

  //scale and translate to [0 , 1]
  for (unsigned int dim = 0; dim<3; dim++) {
    for( size_t ii=0; ii<v.size(); ii++) {
      mn[dim]= std::min(v[ii][dim],mn[dim]);
      mx[dim] = std::max(v[ii][dim],mx[dim]);
    }
    scalar translate = -mn[dim];
    for(size_t ii=0; ii<v.size(); ii++) {
      v[ii][dim]=(v[ii][dim]+translate);
    }
  }

  scalar scale = 1/(mx[0]-mn[0]);
  for(unsigned int dim=1; dim<3; dim++) {
    scale=std::min(1/(mx[dim]-mn[dim]),scale);
  }

  for(size_t ii=0; ii<v.size(); ii++) {
    for (unsigned int dim = 0; dim<3; dim++) {
      v[ii][dim]=v[ii][dim]*scale;
    }
  }
}

void TrigMesh::compute_norm()
{
  n.resize(v.size(), Vector3s::Zero());
  for(unsigned int ii=0; ii<t.size(); ii++) {
    Vector3s a = v[t[ii][1]] - v[t[ii][0]];
    Vector3s b = v[t[ii][2]] - v[t[ii][0]];
    a.cross(b);
    b.normalize();
    for(int jj=0; jj<3; jj++) {
      n[t[ii][jj]]+=b;
      if(t[ii][jj]>=(int)n.size() || t[ii][jj]<0){
        std::cout<<ii<<" "<<jj<<" "<<t[ii][jj]<<" normal computation error\n";
      }
    }
  }
  for(unsigned int ii=0; ii<v.size(); ii++) {
    n[ii].normalize();
  }
}
