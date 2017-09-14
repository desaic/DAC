#include "World.hpp"
#include "Contact.hpp"
#include "CCDTrigMesh.hpp"
#include "EigenUtil.hpp"
#include "Element.hpp"
#include "ElementMesh.hpp"
#include "ElementMeshUtil.hpp"
#include "Material.hpp"
#include "PtTrigIntersect.hpp"
#include "Stepper.hpp"
#include "TrigMesh.hpp"
#include "Timer.hpp"
#include "ArrayUtil.hpp"
#include "FileUtil.hpp"
#include "DeformModel.h"
#include <thread>
#include <set>
#include <string>

using namespace std;

template <typename T>
struct PointCloud
{
  struct Point
  {
    T  x, y, z;
  };

  std::vector<Point>  pts;

  // Must return the number of data points
  inline size_t kdtree_get_point_count() const { return pts.size(); }

  // Returns the distance between the vector "p1[0:size-1]" and the data point with index "idx_p2" stored in the class:
  inline T kdtree_distance(const T *p1, const size_t idx_p2, size_t /*size*/) const
  {
    const T d0 = p1[0] - pts[idx_p2].x;
    const T d1 = p1[1] - pts[idx_p2].y;
    const T d2 = p1[2] - pts[idx_p2].z;
    return d0*d0 + d1*d1 + d2*d2;
  }

  // Returns the dim'th component of the idx'th point in the class:
  // Since this is inlined and the "dim" argument is typically an immediate value, the
  //  "if/else's" are actually solved at compile time.
  inline T kdtree_get_pt(const size_t idx, int dim) const
  {
    if (dim == 0) return pts[idx].x;
    else if (dim == 1) return pts[idx].y;
    else return pts[idx].z;
  }

  // Optional bounding-box computation: return false to default to a standard bbox computation loop.
  //   Return true if the BBOX was already computed by the class and returned in "bb" so it can be avoided to redo it again.
  //   Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3 for point clouds)
  template <class BBOX>
  bool kdtree_get_bbox(BBOX& /*bb*/) const { return false; }
};

World::World():
  contact_h(1e-5f),contactModel(0),handle_contact(1),
  flight_h(1e-4f), flight_fmax(1e-6f),
  terminateOnCollision(0),
  activeMesh(0), stepperIdx(0),
  frameCnt(0),loading_t(0),
  stage(SIM_LOADING),stepper(0),
  frameSkip(1), savePos(0)
{
}

World::~World()
{
  for(unsigned int ii=0; ii<em_.size(); ii++){
    delete em_[ii];
  }
  for(unsigned int ii=0; ii<stepperList.size(); ii++){
    delete stepperList[ii];
  }
  for(unsigned int ii=0; ii<trigm.size(); ii++){
    delete trigm[ii];
  }
  for(unsigned int ii=0; ii<materials.size(); ii++){
    delete materials[ii];
  }

}

bool detectCollision(ElementMesh * m,World * world,
                     const std::vector<Eigen::Vector3d> & x0,
                     std::vector<Contact> & contact)
{
  contact.clear();
  for (size_t mi = 0; mi < world->trigm.size(); mi++){
    updateSurfVert(m->x, world->trigm[mi], world->trigm[mi]->vidx);
  }

  bool hasCollision = false;
  world->collideWith = 0;
  int dim = m->dim;
  for(unsigned int ii = 0; ii<m->x.size(); ii++){
    for(int jj = 0; jj<dim; jj++){
      if(m->fixedDof[dim * ii + jj]){
        continue;
      }
      if( (m->x[ii][jj] < m->lb[ii][jj]) && (m->x[ii][jj] - x0[ii][jj]<0)) {
        hasCollision = true;
        Contact c;
        c.m[0] = 0;
        c.m[1] = -1;
        c.v1 = ii;
        c.normal = Eigen::Vector3d(0, 1, 0);
        c.x = m->x[ii];
        contact.push_back(c);
        world->collideWith = 1;
        //std::cout << "collision " << ii <<" "<<jj<<" "<<m->x[ii][jj]<<" "<<m->lb[ii][jj]<< "\n";
      }
    }
  }

  //box collision
  double eps = 1e-4;
  world->boxContactNormal = Eigen::Vector3d::Zero();
  for(unsigned int ii = 0; ii<m->x.size(); ii++){
    bool collide = true;
    
    for(int jj = 0; jj<dim; jj++){
      if( m->x[ii][jj] < world->box.mn[jj] || m->x[ii][jj] > world->box.mx[jj]) {
        collide = false;
      }
    }
    
    if(collide){
      world->boxContactNormal = Eigen::Vector3d(-1, 0, 0);
      double dx = m->x[ii][0] - world->box.mn[0];
      double dx1 = world->box.mx[0] - m->x[ii][0];
      double dy = world->box.mx[1] - m->x[ii][1];
      //std::cout << "dx dy" << dx << " " << dy << "\n";
      if ( (dx < dy && dx>=0 ) || (dx1<dy && dx1>=0) ){
        if (dx1 < dx && dx1 >= 0){
          world->boxContactNormal = Eigen::Vector3d(1, 0, 0);
        }
        else{
          world->boxContactNormal = Eigen::Vector3d(-1, 0, 0);
        }
      }
      else{
        world->boxContactNormal = Eigen::Vector3d(0, 1, 0);
      }

      world->collideWith = 2;
      Contact c;
      c.m[0] = 0;
      c.m[1] = -1;
      c.v1 = ii;
      c.normal = world->boxContactNormal;
      //c.normal = Eigen::Vector3d(-1, 0, 0);
      c.x = m->x[ii];
      contact.push_back(c);
    }
  }

  //box vertex against mesh
  int trigv[2][3] = { { 0, 1, 2 }, { 0, 2, 3 } };
  for (unsigned int i = 0; i < world->box.edgeVerts.size(); i++){
    Eigen::Vector3d bv = world->box.edgeVerts[i];
    for (int j = 0; j < m->e.size(); j++){
      Element * e = m->e[j];
      
      //loop over hex faces
      for (int f = 0; f < 6; f++){
        //loop over triangles
        for (int t = 0; t < 2; t++){
          Eigen::Vector3d v[3], v0[3];
          for (int vi = 0; vi < 3; vi++){
            int vidx = e->at(HexFaces[f][trigv[t][vi]]);
            v[vi] = m->x[vidx];
            v0[vi] = x0[vidx];
          }

          InterResult result = PtTrigIntersect(bv, bv, v0, v);
          if (result.intersect){
            world->collideWith = 2;
            for (int vi = 0; vi < 4; vi++){
              int vidx = e->at(HexFaces[f][vi]);
              Contact c;
              c.m[0] = 0;
              c.m[1] = -1;
              c.v1 = vidx;
              c.normal = -result.n;
              c.x = (1 - result.t) * x0[vidx] + result.t * m->x[vidx];
              contact.push_back(c);
            }
            break;
          }
        }
      }
    }
  }

  
  if (world->trigm.size() > 0){
    //self
    std::cout << "World.cpp detect self collision.\n";
    world->ccd.update((double*)(world->trigm[0]->v.data()), world->trigm[0]->v.size());
    for (size_t i = 0; i < world->ccd.mdl->vflist.size(); i++){
      int fid = world->ccd.mdl->vflist[i].f;
      int vid = world->ccd.mdl->vflist[i].v;
      TrigMesh * tm = world->trigm[0];
      ElementMesh * em = world->em_[0];
      std::cout << fid << " " << vid <<"\n";
      Contact c;
      c.m[0] = 0;
      c.m[1] = 0;
      //vertex index in fem mesh.
      int ev = tm->vidx[vid];
      c.v1 = ev;
      Eigen::Vector3d trig[3];
      Eigen::Vector3d trigNormal;
      for (int l = 0; l < 3; l++){
        int vi = tm->vidx[tm->t[fid][l]];
        trig[l] = x0[vi];
      }
      trigNormal = (trig[1] - trig[0]).cross(trig[2] - trig[0]);
      trigNormal.normalize();
      c.normal = trigNormal;
      c.x = x0[ev];

      //quad vertex indices in the finite element mesh.
      for (int l = 0; l < 4; l++){
        int vidx = em->e[tm->ei[fid]]->at(HexFaces[tm->fi[fid]][l]);
        c.v2[l] = vidx;
      }

      //corner of quad.
      Eigen::Vector3d x0 = em->x[c.v2[0]];
      Eigen::Vector3d x1 = em->x[c.v2[1]];
      Eigen::Vector3d x2 = em->x[c.v2[2]];

      Eigen::Vector3d a1 = (x1 - x0);
      double l1 = a1.norm();
      a1 /= l1;
      Eigen::Vector3d a2 = (x2 - x0);
      double l2 = a2.norm();
      Eigen::Vector3d a3 = a1.cross(a2).normalized();
      a2 = a3.cross(a1);
      Eigen::Vector3d disp = c.x - x0;
      c.alpha[0] = (float)(a1.dot(disp) / l1);
      c.alpha[1] = (float)(a2.dot(disp) / l2);
      contact.push_back(c);
      std::cout << "contact alpha " << c.alpha[0] << " " << c.alpha[1] << "\n";
      std::cout << "contact N " << c.normal[0] << " " << c.normal[1] << " "<<c.normal[2] << "\n";
      contact.push_back(c);
    }
  }

  if (contact.size()>0){
    std::cout << contact.size() << " contacts.\n";
    //world->trigm[0]->save_obj("debug.obj");
    hasCollision = true;
  }
  return hasCollision;
}

void World::resolveCollisionVer0(std::vector<Contact> & contact)
{
  ElementMesh*em = em_[activeMesh];
  int dim = em->dim;
  std::vector<int> fixedContraints = em->fixedDof;
  for (size_t i = 0; i < contact.size(); i++){
    //infinite friction.
    //assuming collision with y=0 plane.
    int vidx = contact[i].v1;
    for(int jj =0; jj<dim; jj++){
      em->fixedDof[vidx*dim + jj] = 1;
    }
  }
  int ret = stepper->stepWrapper();
  em->fixedDof = fixedContraints;
}

void World::resolveCollisionVer1(std::vector<Contact> & contact)
{
  Stepper * stepper = stepperList[stepperIdx];
  stepper->resolveCollision(contact);
}

double infNorm(const std::vector<Eigen::Vector3d> & v)
{
  double n = 0;
  for(unsigned int ii = 0; ii<v.size(); ii++){
    for(int jj = 0; jj<v[ii].size(); jj++){
      n = std::max(n, std::abs(v[ii][jj]));
    }
  }

  return n;
}

double L1Norm(const std::vector<Eigen::Vector3d> & v)
{
  double n = 0;
  for(unsigned int ii = 0; ii<v.size(); ii++){
    for(int jj = 0; jj<v[ii].size(); jj++){
      n += std::abs(v[ii][jj]);
    }
  }

  return n;
}


double dot(std::vector<Eigen::Vector3d> & a, std::vector<Eigen::Vector3d> & b)
{
  double prod = 0;
  for (size_t i = 0; i < a.size(); i++){
    prod += a[i].dot(b[i]);
  }
  return prod;
}

void dropMesh(ElementMesh * em, World * w, float h0)
{
  std::vector<Eigen::Vector3d> x1 = em->x;

  for (int i = 0; i < 500; i++){
    float dt = 1e-3;
    float h = h0 - 0.5 * 9.8 * std::pow(dt * i, 2);
    if (h < 0){
      break;
    }
    for (size_t i = 0; i < em->x.size(); i++){
      em->x[i][1] = x1[i][1] + h;
    }
    w->saveMeshState();
    w->frameCnt++;
  }
  exit(0);
}

void World::loop()
{
  std::vector<Eigen::Vector3d> x0,v0, vu0;
  ElementMesh * em = em_[activeMesh];
  std::vector<Contact> contact;
  stage = SIM_LOADING;
  double totalTime = 0;
  frameCnt = 0;
  totalAng = 0;
  //double Ene = em->getAllEnergy();
  //std::cout << "Energy: " << Ene << "\n";
  Eigen::Vector3d cen = em->centerOfMass();
  Eigen::Vector3d L0 = em->angularMomentum(cen);
  Eigen::Vector3d center = em->centerOfMass();
  height = center[1];
  //Eigen::AngleAxisd R = Eigen::AngleAxisd(0.02, Eigen::Vector3d(0, 0, 1));
  //Eigen::AngleAxisd R = Eigen::AngleAxisd(3.14, Eigen::Vector3d(0, 1, 0));
  //Eigen::AngleAxisd R1 = Eigen::AngleAxisd(0.06, Eigen::Vector3d(0, 0, 1));
  //Eigen::AngleAxisd R = Eigen::AngleAxisd(-2.51, Eigen::Vector3d(0, 0, 1));
  //Eigen::AngleAxisd R1 = Eigen::AngleAxisd(0, Eigen::Vector3d(0, 0, 1));
  //Eigen::AngleAxisd R = Eigen::AngleAxisd(0.4, Eigen::Vector3d(0, 0, 1));
  //Eigen::AngleAxisd R1 = Eigen::AngleAxisd(0.05, Eigen::Vector3d(1, 0, 0));

  for (size_t i = 0; i < em->x.size(); i++){
    //em->x[i] = R1*R*em->x[i];
    //em->x[i][1] += 0.03;
    //em->v[i][1] = -1;
    //em->v[i][1] = -0.767;
    //em->v[i][1] = 0;
    //if (em->x[i][1] > 0.03){
    //  em->v[i][2] += 0.05;
    //}
  }
  //Eigen::Vector3d mx(0, 0, 0), mn(100, 100, 100);
  //BBox(em->x, mn, mx);
  //std::cout <<"mn "<< mn[1] << "\n";

  float eps = 1e-5f;
  //fixed vertices
  //for (unsigned int ii = 0; ii < em->x.size(); ii++){
  //  if (em->x[ii][1] < mn[1] + eps){
  //    //&& em->X[ii][0] < min[0] + 0.036){
  //    for (int jj = 0; jj < 3; jj++){
  //      em->fixedDof[3 * ii + jj] = 1;
  //    }
  //  }
  //}

  //for (size_t i = 0; i < em->x.size(); i++){
  //  em->x[i][1] = em->x[i][1] - mn[1] + 1e-5f;
  //}

  for (size_t i = 0; i < em->v.size(); i++){
    //em->v[i][0] = 1;
  }

  if (trigm.size()>0 && handle_contact){
    ccd.init((double*)(trigm[0]->v.data()), trigm[0]->v.size(), (int*)(trigm[0]->t.data()), trigm[0]->t.size());
  }

  //make feet soft
  float feet_height = 0.003;
  //for (size_t i = 0; i < em->e.size(); i++){
  //  if (em->X[em->e[i]->at(2)][1] < feet_height - eps){
  //    em->me[i] = 1;
  //  }
  //}

  //dropMesh(em, this, 0.03f);
  //saveStiffness(em);
  //saveModalMeshes(em);
  //exit(0);
  //saveMeshState();
  for(stepperIdx = 0; stepperIdx<stepperList.size(); stepperIdx++){
    stepper = stepperList[stepperIdx];
    stepper->handle_contact = handle_contact;
    //for loading of fixed objects such calibration rigs.
    //if (stepperIdx == 0){
    //  stepper->handle_contact = false;
    //}

    stepper->solvers = &solvers;
    stepper->m=em;
    stepper->meshes = em_;
    if (stepper->simType != Stepper::SIM_RIGID){
      stepper->solver = solvers[activeMesh];
    }    
    stepper->init();
    
    if(stepperIdx>0){
      int nAng = (int)stepperList[stepperIdx-1]->angSeq.size();
      if (nAng>0){
        float ang0 = stepperList[stepperIdx - 1]->angSeq[nAng - 1];
//        std::cout<<ang0<<" angle 0\n";
        stepper->angSeq.push_back(ang0);
      }
    }
   
    stepper->h = contact_h;
    int dim = em->dim;
    //if (stepperIdx == 2){
    //  for (size_t i = 0; i < em->x.size(); i++){
    //    if (em->X[i][0]<0.01)
    //    em->lb[i][1] = -0.001;
    //  }
    //}
    for(int iter = 0; iter<stepper->nSteps; iter++){
      frameCnt++;
      std::cout<<"frame "<<frameCnt<<"\n";
      //predict
      x0 = em->x;
      v0 = em->v;

      Timer timer;
      timer.startWall();
      int ret = stepper->stepWrapper();
      timer.endWall();
      std::cout << "stepper time " << timer.getSecondsWall() << "\n";

      cen = em->centerOfMass();
      Eigen::Vector3d L = em->angularMomentum(cen);
      std::cout << "L " << L[2] << "\n";
      bool correctL = false;
      if (correctL &&  (std::abs(L0[2])*0.999 - std::abs(L[2]) > 0 ) ){
        Eigen::Vector3d p = em->linearMomentum();
        double M = 0;
        for (unsigned int i = 0; i < em->mass.size(); i++){
          M += em->mass[i];
        }

        Eigen::Vector3d vl = (1.0 / M) * p;

        for (size_t i = 0; i < em->v.size(); i++){
          Eigen::Vector3d x = em->x[i] - cen;
          if (x.squaredNorm() < 1e-20){
            continue; 
          }
          x.normalize();
          Eigen::Vector3d v = em->v[i];
          //remove linear part.
          v -= vl;
          Eigen::Vector3d L_i = x.cross(v);
          Eigen::Vector3d angularComp = L_i.cross(x);
          v = v + ((L0[2] - L[2]) / L[2]) * angularComp;
          //add back linear part.
          v += vl;
          em->v[i] = v;
        }
        L = em->angularMomentum(cen);
        std::cout << "L after " << L[2] << "\n";
        std::cout << ret << "\n";
        p = em->linearMomentum();
        std::cout << "p after " << p[1] << "\n";
      }

      if(ret>0){
        if( stepper->simType != Stepper::SIM_DYNAMIC){
          if(handle_contact){
            std::fill(em->fixedDof.begin(), em->fixedDof.end(), 0);
          }
          std::fill(em->fe.begin(), em->fe.end(), Eigen::Vector3d::Zero());
          if(stage==SIM_LOADING){
            stage = SIM_LAUNCHING;
          }
        }
        saveMeshState();
        break;
      }
      //saveMeshState();
      //frameCnt++;
      //if( (em->specialVert >= 0)
      //  && (stage == SIM_LOADING)
      //  && (stepper->simType == Stepper::SIM_DYNAMIC)
      //  && (em->x[em->specialVert][1] <= em->lb[em->specialVert][1]) ){
      //  stage = SIM_LAUNCHING;
      //  std::cout<<"Loading collision at time "<<totalTime<<".\n";
      //  for(unsigned int ii =0; ii<em->x.size(); ii++){
      //    em->fe[ii] = Eigen::Vector3d::Zero();
      //    em->v[ii]  = Eigen::Vector3d::Zero();
      //  }
      //  em->lb[em->specialVert][1] = 0;
      //}

      bool hasCollision = false;
      if(handle_contact
         && stepper->simType != Stepper::SIM_STATIC){
        contact.clear();
        hasCollision = detectCollision(em, this, x0, contact);
      }
      if(!hasCollision){
        totalTime += stepper->h;
        em->vu.clear();
      }
//      std::cout<<"Collision "<<hasCollision<<" "<<stage<<"\n";
      //add check min>some epsilon.
      if (!hasCollision && stage == SIM_LAUNCHING && totalTime > 0.02){// && totalTime > 0.015){
        stage = SIM_FLIGHT;
        stepper->h = flight_h;
        std::cout<<"Launch time "<<totalTime<<"\n";
      }
      if(stage == SIM_FLIGHT && stepperIdx<stepperList.size()-1
         && stepper->simType == Stepper::SIM_DYNAMIC){
          break;
      }
      if(hasCollision){
        stepper->h = contact_h;
      }else if(stage == SIM_FLIGHT){
        stepper->h = flight_h;
      }

      timer.startWall();
      //collision resolution
      if(hasCollision){
        //save contact
        if (contact.size() > 0){
          
          std::string prefix = simDirName + "/" + "c";
          std::string filename = sequenceFilename(frameCnt, 0, prefix.c_str());
          std::ofstream out(filename);
          if (!out.good()){
            //std::cerr<<"cannot open output file "<<filename<<"\n";
          }
          else{
            saveContact(out, contact);
          }
          out.close();
        }
        if( stepper->simType != Stepper::SIM_RIGID){
          em->x = x0;
          em->v = v0;
          if(contactModel==0){
            resolveCollisionVer0(contact);
          }else{
            std::cout << "resolve "<<contact.size()<<" collisions\n";
            resolveCollisionVer1(contact);
          }
        }
        totalTime += stepper->h;
        Eigen::Vector3d p = em->linearMomentum();
        if( ((terminateOnCollision && stage==SIM_FLIGHT)
             || stepper->simType == Stepper::SIM_RIGID)
            && ((p[1]<0 && collideWith==1) || collideWith==2)) {
          em->x = x0;
          em->v = v0;
          std::cout<<"time: "<<totalTime<<" p_y: "<<p[1] <<"\n";
          break;
        }
      }
      timer.endWall();
      std::cout << "Collision time " << timer.getSecondsWall() << "\n";
      if (savePos &&  frameCnt % frameSkip == 0){
        saveMeshState();
      }
      center = em->centerOfMass();
      if (center[1]>height){
        height = center[1];
      }
    }
    std::cout << "Exit stepper loop "<<stepperIdx<<"\n";
    saveMeshState();
    double ang;
    ang = em->getTopAng();
    //Ene = em->getAllEnergy();
    //std::cout << "Energy: " << Ene << "\n";
    std::cout<<"current ang "<<ang<<"\n";
    std::cout << "height " << height << "\n";
    stepper->angSeq.push_back((float)ang);
    ang= stepper->totalAng();
    totalAng += ang;
//    std::cout<<ang<<" angle\n";
    bool savefile = true;
    stepper->finalize(savefile);
  }
  center = em->centerOfMass();
  //kinetic.
  double Ek = em->getKineticEnergy(true);
  //gravity potential
  
  Eigen::Vector3d L = em->angularMomentum(center);
  int NFlip = (int)(- totalAng / (2 * M_PI) );
  double landAng = -totalAng - NFlip * 2 * M_PI;
  if (landAng > M_PI){
    landAng -= 2 * M_PI;
  }

  if (em_.size() > 0){
    Eigen::Vector3d mn, mx;
    BBox(em_[0]->x, mn, mx);
    std::cout << "final bounding box\n";
    std::cout << mn[0] << " " << mn[1] << " " << mn[2] << "\n";
    std::cout << mx[0] << " " << mx[1] << " " << mx[2] << "\n";
  }

  std::cout<<"Total angle "<<totalAng<<" "<<landAng<<"\n";
  std::cout << "collide with " << collideWith << "\n";
  std::cout << "normal " << boxContactNormal[0] <<" "<<boxContactNormal[1]<< "\n";
  std::cout << "Ek " << Ek << "\n";
}

void World::saveMeshState()
{
  std::string prefix = simDirName + "/" + std::to_string(activeMesh) + "m";
  std::string filename = sequenceFilename(frameCnt, NULL, prefix.c_str());
  std::ofstream out(filename);
  if(!out.good()){
//    std::cerr<<"cannot open output file "<<filename<<"\n";
    return;
  }else{
    save_arr(em_[activeMesh]->x, out);
    save_arr(em_[activeMesh]->v, out);
    save_arr(em_[activeMesh]->fe, out);
    save_arr(em_[activeMesh]->fc, out);
  }
  out.close();
}

void runWorld(World * world)
{
  world->loop();
}