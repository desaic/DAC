#include "WorldC.hpp"
#include "ConfigFile.hpp"
#include "Contact.hpp"
#include "Element.hpp"
#include "ElementMesh.hpp"
#include "ElementMeshUtil.hpp"
#include "FileUtil.hpp"
#include "LinSolver.hpp"
#include "LinPardiso.hpp"
#include "ConfigUtil.hpp"
#include "MaterialQuad.hpp"
#include "StrainEne.hpp"
#include "Quadrature.hpp"
#include "ArrayUtil.hpp"
#include "StepperNewmark.hpp"
#include "StepperStatic.hpp"
#include "StepperRigid.hpp"
#include "Stepper.hpp"
#include "StretchMesh.hpp"
#include "PtTrigIntersect.hpp"

WorldC::WorldC() :
staticSim(0),
cubeLocation(0,0,0),
contact_h(1e-5f), contactModel(0), handle_contact(1),
hasInitPos(false),
flight_h(1e-4f),
launchAreaSize(1.0f),
terminateOnCollision(0),
activeMesh(0), stepperIdx(0),
frameCnt(0),
stepper(0),
frameSkip(1), savePos(0)
{
  staticSim = new StepperStatic();
}

WorldC::~WorldC()
{
  for (unsigned int ii = 0; ii<em_.size(); ii++){
    delete em_[ii];
  }
  for (unsigned int ii = 0; ii<stepperList.size(); ii++){
    delete stepperList[ii];
  }
  for (unsigned int ii = 0; ii<trigm.size(); ii++){
    delete trigm[ii];
  }
  for (unsigned int ii = 0; ii<materials.size(); ii++){
    delete materials[ii];
  }
}

void WorldC::resolveCollision(std::vector<Contact> & contact)
{
  Stepper * stepper = stepperList[stepperIdx];
  stepper->meshes = em_;
  stepper->resolveCollision(contact);
}

///@param dist. Distance in meters from top left front corner, not mn or mx.
///to accomendate parametric shape changes.
///@return index of element whose vertex 2 is closest to the fixture point.
///@param w float[4] bilinear weights on the quad.
int findFixtureFrame(ElementMesh * em, Eigen::Vector3d dist, float * w)
{
  //find element whose e[2] is closest to the desired location specified by dist.
  //and its x z coordinates are less than coordinate of the point
  int ev = 2;
  Eigen::Vector3d mn, mx;
  BBox(em->x, mn, mx);
  double minDist = 0;
  int minEle = -1;
  dist[0] = mn[0] + dist[0];
  dist[1] = mx[1] - dist[1];
  dist[2] = mn[2] + dist[2];
  for (size_t i = 0; i < em->e.size(); i++){
    Eigen::Vector3d v = em->x[em->e[i]->at(ev)];
    if (v[0]>dist[0] || v[2] > dist[2]){
      continue;
    }
    double d = (v - dist).norm();
    if (minEle < 0 || d < minDist){
      minDist = d;
      minEle = (int)i;
    }
  }
  float alpha[2];
  Eigen::Vector3d v = em->x[em->e[minEle]->at(ev)];
  double dx = em->x[em->e[minEle]->at(6)][0] - v[0];
  double dy = em->x[em->e[minEle]->at(3)][2] - v[2];
  alpha[0] = (float)((dist[0] - v[0]) / dx);
  alpha[1] = (float)((dist[2] - v[2]) / dy);
  bilinearWeights(alpha, w);
  return minEle;
}

void WorldC::solveStatic()
{
  int nStaticSteps = 500;
  staticSim->m = em_[0];
  staticSim->handle_contact = handle_contact;
  staticSim->solver = solvers[0];
  staticSim->init();
  staticSim->fixed0 = em_[0]->fixedDof;
  std::cout << "static sim\n";
  for (int iter = 0; iter < nStaticSteps; iter++){
    int ret = staticSim->oneStep();
    if (savePos){
      std::cout << iter << ": " << ret << "\n";
      saveMeshState();
    }
    frameCnt++;
    if (ret != 0){
      break;
    }
  }
  em_[0]->fixedDof = staticSim->fixed0;
}

void WorldC::initSolvers()
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

bool intersectBox(const Box & box, const std::vector<Eigen::Vector3d> & x)
{
  bool intersect = false;
  int dim = 3;
  for (size_t i = 0; i < x.size(); i++){
    bool inside = true;
    for (int d = 0; d < dim; d++){
      if (x[i][d]<box.mn[d] || x[i][d]>box.mx[d]){
        inside = false;
        break;
      }
    }
    if (inside){
      intersect = true;
      break;
    }
  }
  return intersect;
}

void WorldC::loop()
{
  std::vector< std::vector<Eigen::Vector3d> > x0(em_.size()), v0(em_.size());
  int dim = em_[0]->dim; 
  double totalTime = 0;
  stage = LAUNCH;
  if (!hasInitPos){
    solveStatic();
    std::cout << "Finish static sim.\n";
    if (em_.size() > 1 && cubeFixture.eidx >= 0){
      int ev[4] = { 2, 6, 3, 7 };
      Eigen::Matrix3d frame;
      Eigen::Vector3d V[4];
      Eigen::Vector3d Origin = Eigen::Vector3d::Zero();
      for (int i = 0; i < 4; i++){
        V[i] = em_[0]->x[em_[0]->e[cubeFixture.eidx]->at(ev[i])];
        Origin += cubeFixture.bw[i] * V[i];
      }
      //move upwards a little
      double offset = 0.02*(V[1][0] - V[0][0]);
      Origin[1] += offset;
      Eigen::Vector3d x = (V[1] - V[0]).normalized();
      Eigen::Vector3d z = (V[2] - V[0]).normalized();
      Eigen::Vector3d y = z.cross(x);
      frame.col(0) = x;
      frame.col(1) = y;
      frame.col(2) = z;
      applyFrame(em_[1], Origin, frame);
    }
  }
  ElementMesh *  em = em_[0];
  em->fe = std::vector<Eigen::Vector3d>(em->x.size(), Eigen::Vector3d::Zero());
  stepperIdx = 0;
  
  //dynamic launching 

  StepperNewmark * stepper = (StepperNewmark *)stepperList[stepperIdx];
  stepper->handle_contact = handle_contact;
  stepper->solvers = &solvers;
  stepper->init();

  //std::ofstream outf("L.txt");
  stepper->nSteps = 0;
  for (int iter = 0; iter < stepper->nSteps; iter++){
    stepper->h = contact_h;
    //predict
    activeMesh = 0;

    //free flight only simulates the box.
    if (stage == FLIGHT){
      activeMesh = 1;
    }
    for (; activeMesh < em_.size(); activeMesh++){
      ElementMesh *  em = em_[activeMesh];
      stepper->solver = solvers[activeMesh];

      std::vector<bool> collision(em->X.size() * dim);
      stepper->m = em;
      x0[activeMesh] = em->x;
      v0[activeMesh] = em->v;
      int ret = stepper->stepWrapper();
    }
    Eigen::Vector3d centerOfMass = em_[0]->centerOfMass();
    Eigen::Vector3d L = em_[0]->angularMomentum(centerOfMass);
    //outf << L[0] << " " << L[1] << " " << L[2] << "\n";
    if (savePos){
      saveMeshState();
    }
    bool hasCollision = false;
    //detect collision only during launching.
    if (handle_contact && stage == LAUNCH){
      hasCollision = detectCollision(em_, trigm, x0, contact);
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

    }
    if (hasCollision){
      stepper->h = contact_h;
      //restore position
      for (size_t i = 0; i < em_.size(); i++){
        em_[i]->x = x0[i];
        em_[i]->v = v0[i];
      }
      if (stepper->simType != Stepper::SIM_RIGID){
        resolveCollision(contact);
        //if (savePos){
        //  saveMeshState();
        //  frameCnt++;
        //}
#ifdef _DEBUG
        //check if new collisions are created
        std::vector<Contact> debugContact;
        bool hasNewCollision = detectCollision(em_, trigm, x0, debugContact);
        if (hasNewCollision){
          std::cout << "Contact resolution creates new collision " << frameCnt << "\n";
        }
#endif
      }
    }
    totalTime += stepper->h;
    frameCnt++;
    if (stage == LAUNCH && em_.size() > 1){
      //make the box slightly larger.
      Eigen::Vector3d momentum = em_[0]->linearMomentum();
      //Box launchRegion = catapultBox;
      //launchRegion.mx *= launchAreaSize;
      //bool insideBox = intersectBox(launchRegion, em_[1]->x);
      //if (!insideBox){
      //if y momentum of the catapult is less than 0, the projectile is considered launched.
      if (momentum[1]<=0){
        std::cout << "Flight " << frameCnt << "\n";
        stage = FLIGHT;
        break;
      }
    }
  }
  //continue a sim.
  //frameCnt = 4045;
  //rigid body projectile free flight
  stage = FLIGHT;
  if (em_.size()>1 && stepperList.size()>1){
    stepperIdx = 1;
    StepperRigid * rigid = (StepperRigid *)stepperList[stepperIdx];
    activeMesh = 1;
    rigid->m = em_[activeMesh];
    rigid->init();
    rigid->h = flight_h;
    ElementMesh *  em = em_[activeMesh];
    for (int iter = 0; iter < rigid->nSteps; iter++){
      x0[activeMesh] = em->x;
      int ret = rigid->stepWrapper();
      if (savePos){
        saveMeshState();
      }
      float worldWall = catapultBox.mx[0] + wallDist;
      bool hit = hitWall(worldWall, em->x, 0, 1);
      bool hitFloor = hitWall(-0.003, em->x, 1, -1);
      if (hit || hitFloor){
        break;
      }
      frameCnt++;
    }
  }

  //outf.close();

  if (em_.size()>1){
    Eigen::Vector3d c = em_[1]->centerOfMass();
    Eigen::Vector3d p = em_[1]->linearMomentum();
    c[0] -= (catapultBox.mx[0] + wallDist);
    c[1] -= targetHeight;
    double w_p = 1e-2;
    double obj = c.squaredNorm() + w_p * p[0];
    std::cout << "hit " << c[0] << " " << c[1] << " " << c[2] << "\n";
    std::cout << obj << " obj\n";
  }
}

void WorldC::resetSim()
{
  frameCnt = 0;
  if (!hasInitPos){
    for (size_t i = 0; i < em_.size(); i++){
      em_[i]->x = em_[i]->X;
      std::fill(em_[i]->v.begin(), em_[i]->v.end(), Eigen::Vector3d::Zero());
      //reset initial external forces
      if (fe[i].size() == em_[i]->fe.size()){
        em_[i]->fe = fe[i];
      }
    }
  }
  
  if (em_.size()>0){
    BBox(em_[0]->X, catapultBox.mn, catapultBox.mx);
  }
  if (em_.size()>1){
    cubeFixture.eidx = findFixtureFrame(em_[0], cubeLocation, cubeFixture.bw);
  }

  ///recompute element Jacobian since rest shape may have changed.
  for (size_t i = 0; i < em_.size(); i++){
    ElementMesh * em = em_[i];
    em->initElements(quadrature);
  }
}

int loadMaterials(WorldC * world, const ConfigFile & conf)
{
  int dim = 3;
  conf.getInt("dim", dim);

  std::string materialFile = "../../models/material2.txt";
  if (conf.hasOpt("materialfile")){
    materialFile = conf.dir + conf.getString("materialfile");
  }

  FileUtilIn matIn(materialFile);
  if (!matIn.good()){
    return -1;
  }
  std::vector<double> densities;
  std::vector<StrainEne*> ene = loadMaterials(matIn.in, densities);
  matIn.close();
  const Quadrature* quadrature = &(Quadrature::Gauss2);
  if (dim == 2){
    quadrature = &(Quadrature::Gauss2_2D);
  }
  world->materials.resize(ene.size());
  for (unsigned int ii = 0; ii<world->materials.size(); ii++){
    MaterialQuad * matq = new MaterialQuad();
    if (dim == 2){
      matq->q = quadrature;
      matq->e.resize(matq->q->x.size());
    }
    for (unsigned int jj = 0; jj < matq->e.size(); jj++){
      matq->e[jj] = ene[ii];
    }
    world->materials[ii] = matq;
    world->materials[ii]->density = densities[ii];
  }
  world->quadrature = quadrature;
  return 0;
}

void loadForce(WorldC* world, const ConfigFile & conf)
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
  float forcedist = 0.01; //m
  conf.getFloat("forcedist", forcedist);
  std::vector<int> forceVertices;
  Eigen::Vector3d max(0, 0, 0), min(100, 100, 100);
  BBox(em->X, min, max);
  float eps = 1e-5;
  for (unsigned int ii = 0; ii<em->e.size(); ii++){
    for (int jj = 0; jj<em->e[ii]->nV(); jj++){
      int vi = em->e[ii]->at(jj);
      Eigen::Vector3d x = em->X[vi];
      if (max[1] - x[1]      < eps
        && (x[0] - min[0]) < forcedist){
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
}

void loadConstraints(WorldC* world, const ConfigFile & conf)
{
  bool fixbottom = true;
  conf.getBool("fixbottom", fixbottom);
  float eps = 1e-5;
  if (world->em_.size() == 0){
    return;
  }
  ElementMesh * em = world->em_[0];
  int dim = em->dim;
  Eigen::Vector3d max(0, 0, 0), min(100, 100, 100);
  BBox(em->X, min, max);
  //fixed vertices of the first mesh assumed to be the catapult
  if (fixbottom ){
    for (unsigned int ii = 0; ii<em->x.size(); ii++){
      if (em->X[ii][1] < min[1] + eps){
        for (int jj = 0; jj<dim; jj++){
          em->fixedDof[dim * ii + jj] = 1;
        }
      }
    }
  }
}

int loadCScene(WorldC * world, const ConfigFile & conf)
{
  int dim = 3;
  conf.getInt("dim", dim);
  
  loadMaterials(world, conf);
 
  
  float meshscale = 1;
  conf.getFloat("meshscale", meshscale);
  
  std::vector<std::string > meshfiles = conf.getStringVector("meshfile");
  std::vector<std::string> posfiles = conf.getStringVector("initpos");
  if (meshfiles.size()== 0){
    std::cout << "Need a mesh file\n";
    return -1;
  }
    
  for (size_t i = 0; i < meshfiles.size(); i++){
    std::string meshfile = conf.dir + meshfiles[i];
    FileUtilIn in(meshfile);
    if (!in.good()){
      continue;
    }
    ElementMesh * em = new ElementMesh();
    em->dim = dim;
    em->load(in.in, meshscale);
    in.close();
    for (size_t j = 0; j < world->materials.size(); j++){
      em->addMaterial(world->materials[j]);
    }
    
    conf.getFloat("damp_alpha", em->dampAlpha);
    conf.getFloat("damp_beta", em->dampBeta);
    em->check();
    world->em_.push_back(em);
    TrigMesh * tm = new TrigMesh();
    hexToTrigMesh(em, tm);
    world->trigm.push_back(tm);
  }

  //put the first mesh the catapult at origin
  if (world->em_.size() > 0){
    placeAtOrigin(world->em_[0]);
  }
  loadForce(world, conf);
  loadConstraints(world, conf);

  //can contain duplicates
  conf.getFloat("dt", world->flight_h);
  conf.getFloat("dtcontact", world->contact_h);
  conf.getBool("handle_contact", world->handle_contact);
  conf.getInt("contactmodel", world->contactModel);

  for (size_t i = 0; i < world->em_.size(); i++){
    ElementMesh * em = world->em_[i];
    em->initElements(world->quadrature);
  }
  world->initSolvers();

  std::vector<Stepper*> stepperList;
  std::vector<std::string> stepperTypes = conf.getStringVector("stepper");
  int frictionModel = 0;
  conf.getInt("frictionModel", frictionModel);
  for (unsigned int ii = 0; ii<stepperTypes.size(); ii++){
    Stepper * stepper = 0;
    if (stepperTypes[ii] == "newmark"){
      StepperNewmark * sn = new StepperNewmark();
      sn->frictionModel = (Stepper::FrictionModel)frictionModel;
      stepper = sn;
    }
    else if (stepperTypes[ii] == "rigid"){
      StepperRigid * st = new StepperRigid();
      stepper = st;
    }
    if (stepper != 0){
      stepperList.push_back(stepper);
      conf.getBool("addgravity", stepper->addGravity);
      stepper->logFileName = conf.getString("logfile");
      conf.getInt("nSteps", stepper->nSteps);
    }
  }

  if (stepperList.size() == 0){
    std::cout << "Error in config: need a time stepper\n";
    return -1;
  }

  std::vector<float> cubev = conf.getFloatVector("cubeLocation");
  if (cubev.size() != 3){
    std::cout << "Need to specify cube location.\n";
    world->cubeLocation = Eigen::Vector3d::Zero();
  }
  else{
    for (size_t i = 0; i < cubev.size(); i++){
      world->cubeLocation[i] = cubev[i];
    }
  }

  StretchOpts opts;
  if (conf.hasOpt("beamwidth")){
    conf.getFloat("beamwidth", opts.beamwidth);
    conf.getFloat("leftmargin", opts.leftmargin);
    opts.offset = Eigen::Vector3d(0.0, 0.0, 0);
    std::vector<float> stretch = conf.getFloatVector("stretch");
    if (stretch.size() == 3){
      opts.offset = Eigen::Vector3d(stretch[0], stretch[1], stretch[2]);
    }
    if (world->em_[0]->dim == 2){
      segmentJumper(*(world->em_[0]), opts);
      stretchJumper(*(world->em_[0]), opts);
    }
    else{
      segmentJumper3d(*(world->em_[0]), opts);
      stretchJumper3d(*(world->em_[0]), opts);
    }
  }

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


  conf.getBool("terminateOnCollision", world->terminateOnCollision);
  conf.getBool("savepos", world->savePos);
  world->wallDist = 1;
  conf.getFloat("wallDist", world->wallDist);
  world->simDirName = conf.dir + "/" + conf.getString("resultdir");
  conf.getInt("frameskip", world->frameSkip);
  conf.getFloat("targetHeight", world->targetHeight);
  conf.getFloat("launchAreaSize", world->launchAreaSize);
  world->stepperList = stepperList;

  return 0;
}

void WorldC::saveMeshState()
{
  int mi = 0;
  if (stage == FLIGHT){
    mi = 1;
  }
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
