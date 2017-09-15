#include "WorldStair.hpp"
#include "ConfigFile.hpp"
#include "Contact.hpp"
#include "Element.hpp"
#include "ElementHex.hpp"
#include "ElementMesh.hpp"
#include "ElementMeshUtil.hpp"
#include "ElementRegGrid.hpp"
#include "FileUtil.hpp"
#include "LinSolver.hpp"
#include "LinPardiso.hpp"
#include "ConfigUtil.hpp"
#include "MaterialQuad.hpp"
#include "StrainEne.hpp"
#include "Quadrature.hpp"
#include "ArrayUtil.hpp"
#include "StepperNewmark.hpp"
#include "Stepper.hpp"
#include "PtTrigIntersect.hpp"

void loadHeadVox(ElementRegGrid * head);

WorldStair::WorldStair() :
dt(1e-5f),
forceTime(1e-2f),
dtheta(-50.f),
stairHeight(0.05f), stairHeightStart(0.24f),
stairWidth(0.05f),
wallDist(1.0f),
activeMesh(0),
frameCnt(0),
frameSkip(1),
savePos(0)
{}

WorldStair::~WorldStair()
{
  for (unsigned int ii = 0; ii<em_.size(); ii++){
    delete em_[ii];
  }
  for (size_t i = 0; i < stepperList.size(); i++){
    if (stepperList[i] != 0){
      delete stepperList[i];
    }
  }
  for (unsigned int ii = 0; ii<materials.size(); ii++){
    delete materials[ii];
  }
}

int WorldStair::resolveCollision(std::vector<Contact> & contact)
{
  stepperList[0]->meshes = em_;
  StepperNewmark * sn = (StepperNewmark*)stepperList[0];
  sn->frictionModel = StepperNewmark::FINITE;
  return sn->resolveCollision3D_old(contact);
}

void WorldStair::initSolvers()
{
  for (size_t i = 0; i < em_.size(); i++){
    ElementMesh * em = em_[i];
    if (solvers.size() <= i){
      LinPardiso<double>  * solver = new LinPardiso<double>();
      solver->init();
      Eigen::SparseMatrix<double> K = em->getStiffnessSparse();
      solver->init(K);
      solvers.push_back(solver);
    }
  }
}

bool detectCollision(ElementMesh * m,
  WorldStair * world,
  const std::vector<Eigen::Vector3d> & x0,
   std::vector<Contact> & contact)
{
  bool hasCollision = false;
  int dim = m->dim;
  double eleSize = m->x[1][2] - m->x[0][2];
  float eps = (float)eleSize;
  for (unsigned int ii = 0; ii<m->x.size(); ii++){
    if (m->fixedDof[dim * ii + 1]){
      continue;
    }
    //int stepN = (int)( m->x[ii][0] / world->stairWidth);
    //stepN = std::max(0, stepN);
    //float lb = world->stairHeightStart - stepN * world->stairHeight;
    
    //if (lb - eps > m->x[ii][1]){
    float lb = world->stairHeightStart - (float)m->x[ii][0] * (world->stairHeight / world->stairWidth);
    lb = std::max(0.0f, lb);
    float norm = world->stairWidth * world->stairWidth + world->stairHeight * world->stairHeight;
    if (norm > 1e-20){
      norm = std::sqrt(norm);
    }
    else{
      norm = 1;
    }
    Eigen::Vector3d normal = Eigen::Vector3d(world->stairHeight / norm, world->stairWidth / norm, 0);
    float x1 = world->stairHeightStart / world->stairHeight * world->stairWidth;
    if (m->x[ii][0] > x1){
      normal = Eigen::Vector3d::UnitY();
    }
    //}
    //float lb = 0;
    if ((m->x[ii][1] < lb) ){//&& (m->x[ii][1] - x0[ii][1] < 0)) {
      hasCollision = true;
      Contact c;
      c.m[0] = 0;
      c.m[1] = -1;
      c.v1 = ii;
      c.normal = normal;

      c.x = m->x[ii];
      contact.push_back(c);
      //std::cout << "collision " << ii <<" "<<jj<<" "<<m->x[ii][jj]<<" "<<m->lb[ii][jj]<< "\n";
    }
  }
  return hasCollision;
}

void rotateVec(std::vector<Eigen::Vector3d> & v, float theta)
{
  Eigen::Matrix3d R;
  R = Eigen::AngleAxisd(theta, Eigen::Vector3d::UnitZ());
  for (size_t i = 0; i < v.size(); i++){
    v[i] = R*v[i];
  }
}

bool WorldStair::finishSim()
{
  double x1 = stairWidth * stairHeightStart / stairHeight;
  for (size_t i = 0; i < em_[0]->x.size(); i++){
    if (em_[0]->x[i][0] > x1){
      return true;
    }
  }
  return false;
}

void WorldStair::loop()
{
  std::vector< std::vector<Eigen::Vector3d> > x0(em_.size()), v0(em_.size());
  int dim = em_[0]->dim; 
  double totalTime = 0;
  ElementMesh *  em = em_[0];
  StepperNewmark * stepper = (StepperNewmark *)stepperList[0];
  stepper->handle_contact = true;
  stepper->solvers = &solvers;
  stepper->init();
  past_finish_line = false;
  for (int iter = 0; iter < stepper->nSteps; iter++){
    stepper->h = dt;
    activeMesh = 0;
    int ret = 0;
    for (; activeMesh < em_.size(); activeMesh++){
      ElementMesh *  em = em_[activeMesh];
      stepper->solver = solvers[activeMesh];
      stepper->newton->maxIters = 1;
      std::vector<bool> collision(em->X.size() * dim);
      stepper->m = em;
      x0[activeMesh] = em->x;
      v0[activeMesh] = em->v;
      int ret = stepper->stepWrapper();
      if (ret == -2){
        past_finish_line = true;
        overall_speed = -100;

        break;
      }
    }
    if (ret == -2){
      break;
    }
    bool hasCollision = false;
    contact.clear();
    hasCollision = detectCollision(em_[0], this, x0[0], contact);
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

    if (hasCollision){
      stepper->h = dt;
      //restore position
      for (size_t i = 0; i < em_.size(); i++){
        em_[i]->x = x0[i];
        em_[i]->v = v0[i];
      }
      int status = resolveCollision(contact);
      if (status == -2){
        past_finish_line = true;
        overall_speed = -100;
        break;
      }
    }
    if (savePos){
      saveMeshState();
    }

    totalTime += stepper->h;
    //std::cout << totalTime << "\n";
    if (finishSim()){
      past_finish_line = true;
      finish_time = totalTime;
      centerEnd = em_[0]->centerOfMass();
      Eigen::Vector3d dispxy = centerEnd - centerStart;
      dispxy[2] = 0;
      overall_speed = dispxy.norm() / totalTime;
      break;
    }
    frameCnt++;
  }

  if (!past_finish_line){
    finish_time = totalTime;
    centerEnd = em_[0]->centerOfMass();
    Eigen::Vector3d dispxy = centerEnd - centerStart;
    dispxy[2] = 0;
    overall_speed = dispxy.norm() / totalTime;
  }
  std::cout << "t v " << totalTime<<" " <<overall_speed;
}

void orientMesh(ElementMesh * em, 
  Eigen::Vector3d x, Eigen::Vector3d y, Eigen::Vector3d z)
{
  Eigen::Vector3d mn, mx;
  BBox(em->x, mn, mx);
  for (size_t i = 0; i < em->x.size(); i++){
    Eigen::Vector3d v = em->x[i] - mn;
    em->x[i] = v[0] * x + v[1] * y + v[2] * z;
  }
}

void WorldStair::resetSim()
{
  frameCnt = 0;
  for (size_t i = 0; i < em_.size(); i++){
    em_[i]->x = em_[i]->X;
    std::fill(em_[i]->v.begin(), em_[i]->v.end(), Eigen::Vector3d::Zero());
  }
  
  //move mesh to starting position
  placeAtOrigin(em_[0]);

  Eigen::Vector3d stairSlope;
  stairSlope << stairWidth, -stairHeight, 0;
  stairSlope.normalize();
  Eigen::Vector3d y (0, 1, 0);
  Eigen::Vector3d z(0, 0, 1);
  y = z.cross(stairSlope);
  y.normalize();

  float tilt = 0.2;
  z += tilt * y;
  z.normalize();
  y = z.cross(stairSlope);
  y.normalize();
  //orientMesh(em_[0], stairSlope, y, z);
  float eps = 1e-5;
  Eigen::Vector3d mn, mx;
  BBox(em_[0]->x, mn, mx);
  float h = stairHeight*(mx[0] - mn[0]) / stairWidth;
  translate(em_[0], Eigen::Vector3d(0, stairHeightStart + 1e-5 + h, 0));
  ///recompute element Jacobian since rest shape may have changed.
  for (size_t i = 0; i < em_.size(); i++){
    ElementMesh * em = em_[i];
    em->initElements(quadrature);
  }
  centerStart = em_[0]->centerOfMass();
}

void loadForce(WorldStair * world, const ConfigFile & conf)
{
  world->fe.resize(world->em_.size());
  if (world->em_.size() == 0){
    return;
  }
  ElementMesh * em = world->em_[0];
  Eigen::Vector3d ff(0, 0, 0);
  if (conf.hasOpt("force")){
    std::vector<float> confForce = conf.getFloatVector("force");
    if (confForce.size() != 3){
      std::cout << "Wrong force dimension in config file\n";
    }
    for (int ii = 0; ii<3; ii++){
      ff[ii] = confForce[ii];
    }
  }
  std::vector<int> forceVertices;
  Eigen::Vector3d max(0, 0, 0), min(100, 100, 100);
  BBox(em->X, min, max);
  float eps = 1e-5f;
  for (unsigned int ii = 0; ii<em->e.size(); ii++){
    for (int jj = 0; jj<em->e[ii]->nV(); jj++){
      int vi = em->e[ii]->at(jj);
      Eigen::Vector3d x = em->X[vi];
      if (max[1] - x[1] < eps){
        forceVertices.push_back(vi);
        //std::cout << vi << "\n";
      }
    }
  }

  float forceScale = 1.0f / forceVertices.size();
  world->fe[0] = em->fe;
  for (unsigned int ii = 0; ii<forceVertices.size(); ii++){
    int vi = forceVertices[ii];
    world->fe[0][vi] += forceScale * ff;
  }
  em->fe = world->fe[0];
}

void makeIshape(ElementRegGrid * em, const std::vector<int> & gridres)
{
  //make middle thinner.
  std::vector<Element*> newEle;
  for (int i = 0; i < gridres[0]; i++){
    for (int j = 0; j < gridres[1]; j++){
      for (int k = 0; k < gridres[2]; k++){
        int eIdx = i * gridres[1] * gridres[2] + j * gridres[2] + k;
        Element * e = em->e[eIdx];
        if (j == 0 || j == gridres[1] - 1 ||
          ((i == gridres[0] / 2 && k == gridres[2] / 2))) {
          newEle.push_back(e);
        }
        else{
          delete e;
        }
      }
    }
  }
  em->e = newEle;
  em->rmEmptyVert();
}

//\attach another grid to the bottom of one grid
void attachBottom(ElementRegGrid * body, ElementRegGrid * part,
  int x0, int z0)
{
  const int nTop = 4;
  int botVerts[nTop] = { 0, 1, 4, 5 };
  int topVerts[nTop] = { 2, 3, 6, 7 };
  //translate part so that the first vertices meet.
  //does not handle rotation.
  int bodye = body->GetEleInd(x0, 0, z0);
  int parte = part->GetEleInd(0, part->ny - 1, 0);
  Eigen::Vector3d bodyv = body->X[body->e[bodye]->at(botVerts[0])];
  Eigen::Vector3d partv = part->X[part->e[parte]->at(topVerts[0])];
  Eigen::Vector3d disp = bodyv - partv;
  for (size_t i = 0; i < part->X.size(); i++){
    part->X[i] += disp;
  }

  int vertOffset = body->X.size();
  //append vertices. Allow duplicates. 
  //duplicates to be removed in the end.
  for (int i = 0; i <= part->nx; i++){
    for (int j = 0; j <= part->ny; j++){
      for (int k = 0; k <= part->nz; k++){
        int vidx = part->GetVertInd(i, j, k);
        body->X.push_back(part->X[vidx]);
      }
    }
  }

  int eOffset = body->e.size();
  //append all elements
  for (size_t i = 0; i < part->e.size(); i++){
    ElementHex * ele = new ElementHex();
    for (int j = 0; j < part->e[i]->nV(); j++){
      (*ele)[j] = part->e[i]->at(j) + vertOffset;
    }
    body->e.push_back(ele);
  }

  for (int i = 0; i < part->nx; i++){
    for (int j = 0; j < part->nz; j++){
      int eidx = part->GetEleInd(i, part->ny-1, j) ;
      Element * ele = body->e[eidx + eOffset];
      int bodyeidx = body->GetEleInd(i+x0, 0, j+z0);
      for (int k = 0; k < nTop; k++){
        int bodyvidx = body->e[bodyeidx]->at(botVerts[k]);
        (*ele)[topVerts[k]] = bodyvidx;
      }
    }
  }
}

void attachTop(ElementRegGrid * body, ElementRegGrid * part,
  int x0, int z0)
{
  const int nTop = 4;
  int botVerts[nTop] = { 0, 1, 4, 5 };
  int topVerts[nTop] = { 2, 3, 6, 7 };
  //translate part so that the first vertices meet.
  //does not handle rotation.
  int bodye = body->GetEleInd(x0, body->ny - 1, z0);
  int parte = part->GetEleInd(0, 0, 0);
  Eigen::Vector3d bodyv = body->X[body->e[bodye]->at(topVerts[0])];
  Eigen::Vector3d partv = part->X[part->e[parte]->at(botVerts[0])];
  Eigen::Vector3d disp = bodyv - partv;
  for (size_t i = 0; i < part->X.size(); i++){
    part->X[i] += disp;
  }

  int vertOffset = body->X.size();
  //append vertices. Allow duplicates. 
  //duplicates to be removed in the end.
  for (int i = 0; i <= part->nx; i++){
    for (int j = 0; j <= part->ny; j++){
      for (int k = 0; k <= part->nz; k++){
        int vidx = part->GetVertInd(i, j, k);
        if (vidx < 0){
          continue;
        }
        body->X.push_back(part->X[vidx]);
      }
    }
  }

  int eOffset = body->e.size();
  //append all elements
  for (size_t i = 0; i < part->e.size(); i++){
    ElementHex * ele = new ElementHex();
    for (int j = 0; j < part->e[i]->nV(); j++){
      (*ele)[j] = part->e[i]->at(j) + vertOffset;
    }
    body->e.push_back(ele);
  }

  for (int i = 0; i < part->nx; i++){
    for (int j = 0; j < part->nz; j++){
      int eidx = part->GetEleInd(i, 0, j);
      if (eidx < 0){
        continue;
      }
      Element * ele = body->e[eidx + eOffset];
      int bodyeidx = body->GetEleInd(i + x0, body->ny - 1, j + z0);
      for (int k = 0; k < nTop; k++){
        int bodyvidx = body->e[bodyeidx]->at(topVerts[k]);
        (*ele)[botVerts[k]] = bodyvidx;
      }
    }
  }

}

void scaleBody(ElementRegGrid * body, WalkerParam * param,
  int * bodyRes, int * legRes)
{
  //scale x
  //reference element index for scaling the left part of body.
  //reference is the point that's not moved.
  int ref_eidx = body->GetEleInd(legRes[0] - 1, 0, 0);
  int ref_vidx = body->e[ref_eidx]->at(4);
  double x0 = body->X[ref_vidx][0];
  double dx = param->legSize[0] / legRes[0];
  for (int i = 0; i < legRes[0]; i++){
    for (int j = 0; j <= body->ny; j++){
      for (int k = 0; k <= body->nz; k++){
        int vidx = body->GetVertInd(i, j, k);
        body->X[vidx][0] = x0 - param->legSize[0] + i * dx;
      }
    }
  }

  ref_eidx = body->GetEleInd(body->nx - legRes[0], 0, 0);
  ref_vidx = body->e[ref_eidx]->at(0);
  x0 = body->X[ref_vidx][0];
  for (int i = 1; i <= legRes[0]; i++){
    for (int j = 0; j <= body->ny; j++){
      for (int k = 0; k <= body->nz; k++){
        int vidx = body->GetVertInd(body->nx - legRes[0] + i, j, k);
        body->X[vidx][0] = x0 + i * dx;
      }
    }
  }

  //scale z direction
  ref_eidx = body->GetEleInd(0, 0, legRes[2] - 1);
  ref_vidx = body->e[ref_eidx]->at(1);
  double z0 = body->X[ref_vidx][2];
  double dz = param->legSize[2] / legRes[2];
  
  for (int i = 0; i <= body->nx; i++){
    for (int j = 0; j <= body->ny; j++){
      for (int k = 0; k < legRes[2]; k++){
        int vidx = body->GetVertInd(i, j, k);
        body->X[vidx][2] = z0 - param->legSize[2] + k * dz;
      }
    }
  }

  ref_eidx = body->GetEleInd(0, 0, body->nz - legRes[2]);
  ref_vidx = body->e[ref_eidx]->at(0);
  z0 = body->X[ref_vidx][2];
  for (int i = 0; i <= body->nx; i++){
    for (int j = 0; j <= body->ny; j++){
      for (int k = 1; k <= legRes[2]; k++){
        int vidx = body->GetVertInd(i, j, body->nz - legRes[2] + k);
        body->X[vidx][2] = z0 + k * dz;
      }
    }
  }

}

void makeWalker(ElementRegGrid * em, WalkerParam * param)
{
  //   y/\  
  //    | 
  //    | 
  //    |______\ x
  //   /
  // z/
  //leg x and z size determines outer layer of body mesh resolution.
  //top body with 4 legs attached below at 4 corners.
  const int dim = 3;
  int legRes[dim];
  //resolution of the central part of the body.
  int bodyRes[dim];
  //size of the body including leg mountings.
  float bodyOuterSize[dim] = {param->bodySize[0] + 2*param->legSize[0],
    param->bodySize[1], param->bodySize[2] + 2 * param->legSize[2] };
  ElementRegGrid leg;
  ElementRegGrid head;
  for (int i = 0; i < dim; i++){
    legRes[i] = std::max(1, (int)(param->legSize[i]/param->dx + 0.5));
    bodyRes[i] = std::max(1, (int)(param->bodySize[i] / param->dx + 0.5));
  }
  //attach a cubic head with equal edge lengths.
  int ny = std::min(bodyRes[0], bodyRes[2]);
  head.resize(bodyRes[0], ny, bodyRes[2]);

  //loadHeadVox(&head);
  Eigen::Vector3d mn, mx;
  BBox(head.X, mn, mx);

  bodyRes[0] += 2 * legRes[0];
  bodyRes[2] += 2 * legRes[2];
  em->resize(bodyRes[0], bodyRes[1], bodyRes[2]);
  std::vector<float> scale(3, 0);

  //scale body.
  int i1 = em->GetVertInd(bodyRes[0] - legRes[0], em->ny, bodyRes[2] - legRes[2]);
  int i2 = em->GetVertInd(legRes[0], 0, legRes[2]);
  Eigen::Vector3d bbox = em->X[i1]- em->X[i2];
  for (int i = 0; i < dim; i++){
    scale[i] = param->bodySize[i] / bbox[i];
  }
  scaleMesh(em, scale);
  
  //scale head
  for (int i = 0; i < dim; i++){
    scale[i] = bbox[i]/(mx[i] - mn[i])*scale[i];
  }
  scale[1] = std::min(scale[0], scale[2]);
  scaleMesh(&head, scale);

  leg.resize(legRes[0], legRes[1], legRes[2]);
  bbox = leg.X[leg.X.size() - 1] - leg.X[0];
  for (int i = 0; i < dim; i++){
    scale[i] = param->legSize[i]/bbox[i];
  }
  scaleMesh(&leg, scale);
  scaleBody(em, param, bodyRes, legRes);
  
 
  attachBottom(em, &leg, 0, 0);
  attachBottom(em, &leg, 0, em->nz - leg.nz);
  attachBottom(em, &leg, em->nx - leg.nx, 0);
  attachBottom(em, &leg, em->nx - leg.nx, em->nz - leg.nz);
  attachTop(em, &head, leg.nx, leg.nz);
  em->rmEmptyVert();
  em->initArrays();

  //loadHeadVox(em);
  for (int ei = em->e.size() - head.e.size(); ei < em->e.size(); ei++){
    em->me[ei] = 1;
  }
  //fix top vertices for modal analysis.
  //for (int i = legRes[0]; i <= em->nx - legRes[0]; i++){
  //  for (int k = legRes[2]; k <= em->nz - legRes[2]; k++){
  //    int vidx = em->GetVertInd(i,em->ny, k);
  //    for (int d = 0; d < dim; d++){
  //      em->fixedDof[dim * vidx + d] = 1;
  //    }
  //  }
  //}

}

void makeBench(ElementRegGrid * em, const std::vector<int> & gridres)
{
  std::vector<Element*> newEle;
  for (int i = 0; i < gridres[0]; i++){
    for (int j = 0; j < gridres[1]; j++){
      for (int k = 0; k < gridres[2]; k++){
        int eIdx = i * gridres[1] * gridres[2] + j * gridres[2] + k;
        Element * e = em->e[eIdx];
        bool corneri = (i == 0 || i == gridres[0] - 1);
        bool cornerk = (k == 0 || k == gridres[2] - 1);
        if ( (j==gridres[1] - 1) || (corneri && cornerk) ){
          newEle.push_back(e);
        }
        else{
          delete e;
        }
      }
    }
  }
  em->e = newEle;
  em->rmEmptyVert();
}

void loadHeadVox(ElementRegGrid * mesh)
{
  std::string filename = "D:/workspace/GitHub/DesignDynamics/models/stair/Minion_head_/minion_vox.txt";
  int nx, ny, nz;
  std::ifstream in(filename);
  in >> nx >> ny >> nz;
  int startIdx = mesh->e.size() - nx * ny * nz;
  
  std::vector<Element*> newE;
  for (int i = 0; i < startIdx; i++){
    newE.push_back(mesh->e[i]);
  }
  int cnt = 0;
  for (int i = 0; i < nx; i++){
    for (int j = 0; j < ny; j++){
      for (int k = 0; k < nz; k++){
        int val = 0;
        in >> val;
        int eidx = startIdx + i * ny * nz + j * nz + k;
        Element * ele = mesh->e[eidx];
        if (val > 0){
          newE.push_back(ele);
          cnt++;
        }
        else{
          delete ele;
        }
      }
    }
  }
  mesh->e = newE;
  mesh->initArrays();
  mesh->rmEmptyVert();
  
}

int loadStairScene(WorldStair * world, const ConfigFile & conf)
{
  int dim = 3;

  world->quadrature = &Quadrature::Gauss2;
  loadMaterials(world->materials, conf, world->quadrature);

  std::vector<int> gridres(2, 3);
  std::vector<std::string> posfiles = conf.getStringVector("initpos");

  ElementRegGrid * em = new ElementRegGrid();
  if (conf.hasOpt("gridres")){
    conf.getIntVector("gridres", gridres);
    em->resize(gridres[0], gridres[1], gridres[2]);
  }
  
  em->dim = dim;

  bool Ishape = false;
  bool bench = false;
  bool walker = true;
  if (Ishape){
    makeIshape(em, gridres);
  }
  else if (bench){
    makeBench(em, gridres);
  }
  else if (walker){
       

    WalkerParam param;
    param.bodySize[0] = 0.023703;
    param.bodySize[1] = 0.0131111;
    param.bodySize[2] = 0.0240401;
    param.legSize[0] = 0.00613194;
    param.legSize[1] = 0.0243221;
    param.legSize[2] = 0.00708386;
    //param.bodySize[0] = 0.02;
    //param.bodySize[1] = 0.01;
    //param.bodySize[2] = 0.02;
    //param.legSize[0] = 0.01;
    //param.legSize[1] = 0.0243221;
    //param.legSize[2] = 0.01;

    //param.bodySize[0] = 0.0279596;
    //param.bodySize[1] = 0.0148753;
    //param.bodySize[2] = 0.0279869;
    //param.legSize[0] = 0.00207999;
    //param.legSize[1] = 0.0228329;
    //param.legSize[2] = 0.00264818;

    param.dx = 0.01;
    makeWalker(em, &param);
    std::cout << "nV " << em->X.size() << "\n";
    //GetVertInd not usable/incorrect after this point.
    //FileUtilOut out("femwalk10x.txt");
    //em->saveMesh(out.out);
    //out.close();
    //TrigMesh tm;
    //hexToTrigMesh(em, &tm);
    //tm.save("femwalk10x.obj");
  }
  //em->initArrays();

  //fix top for testing
  //Eigen::Vector3d mn, mx;
  //BBox(em->X, mn, mx);
  float eleSize = (float)(em->X[1][2] - em->X[0][2]);
  float eps = 0.1f * eleSize;
  //for (size_t i = 0; i < em->X.size(); i++){
  //  if (em->X[i][1] > mx[1] - eps){
  //  //if (em->X[i][1] < mn[1] + eps){
  //    em->fixedDof[3 * i] = 1;
  //    em->fixedDof[3 * i + 1] = 1;
  //    em->fixedDof[3 * i + 2] = 1;
  //  }
  //}

  //assign rigid material to top and bottom;
  //std::fill(em->me.begin(), em->me.end(), 1);
  //for (size_t ei = 0; ei < em->e.size(); ei++){
  //  Eigen::Vector3d x0 = em->X[em->e[ei]->at(0)];
  //  Eigen::Vector3d x2 = em->X[em->e[ei]->at(2)];
  //  if ( (x0[1]<mn[1] + eps) || (x2[1] > mx[1] - eps)) {
  //    em->me[ei] = 0;
  //  }
  //}
  if (conf.hasOpt("meshscale")){
    //meters
    std::vector<float> meshScale(1, 3);
    meshScale = conf.getFloatVector("meshscale");
    scaleMesh(em, meshScale);
  }
  
  std::ofstream out("3d.txt");
  em->saveMesh(out);
  out.close();
  for (size_t j = 0; j < world->materials.size(); j++){
    em->addMaterial(world->materials[j]);
  }
  conf.getFloat("damp_alpha", em->dampAlpha);
  conf.getFloat("damp_beta", em->dampBeta);
  em->check();
  world->em_.push_back(em);
 
  //can contain duplicates
  conf.getFloat("dt", world->dt);

  for (size_t i = 0; i < world->em_.size(); i++){
    ElementMesh * em = world->em_[i];
    em->initElements(world->quadrature);
  }
  
  //saveStiffness(em);
  //exit(0);

  world->initSolvers();
  
  loadForce(world, conf);

  std::vector<Stepper*> stepperList;
  std::vector<std::string> stepperTypes = conf.getStringVector("stepper");
  int frictionModel = 0;
  conf.getInt("frictionModel", frictionModel);
  StepperNewmark * sn = new StepperNewmark();
  sn->frictionModel = (Stepper::FrictionModel)frictionModel;
  sn->nSteps = 0;
  conf.getBool("addgravity", sn->addGravity);
  conf.getInt("nSteps", sn->nSteps);
  world->stepperList.push_back(sn);

  if (posfiles.size() == world->em_.size()){
    for (size_t i = 0; i < posfiles.size(); i++){
      std::string posfilename = conf.dir + "/" + posfiles[i];
      FileUtilIn in(posfilename);
      if (in.good()){
        world->hasInitPos = true;
        loadSimState(world->em_[i] , in.in);
      }
      in.close();
    }
  }

  conf.getBool("savepos", world->savePos);
  world->wallDist = 1;
  conf.getFloat("wallDist", world->wallDist);
  world->simDirName = conf.dir + "/" + conf.getString("resultdir");
  conf.getInt("frameskip", world->frameSkip);
  world->forceTime = 1e-2f;
  conf.getFloat("forceTime", world->forceTime);
  conf.getFloat("stairWidth", world->stairWidth);
  conf.getFloat("stairHeight", world->stairHeight);
  conf.getFloat("stairHeightStart", world->stairHeightStart);
  return 0;
}

void WorldStair::saveMeshState()
{
  int mi = 0;
  for (; mi < (int)em_.size(); mi++){
    std::string prefix = simDirName + "/" + std::to_string(mi) + "m";
    std::string filename = sequenceFilename(frameCnt, 0, prefix.c_str());
    std::ofstream out(filename);
    if (!out.good()){
      //std::cerr<<"cannot open output file "<<filename<<"\n";
      continue;
    }
    else{
      saveSimState(em_[mi], out);
    }
    out.close();
  } 
}
