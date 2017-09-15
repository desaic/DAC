#include "ArrayUtil.hpp"
#include "ConfigUtil.hpp"
#include "FileUtil.hpp"
#include "Element.hpp"
#include "ElementMesh.hpp"
#include "ElementMeshUtil.hpp"
#include "StrainEne.hpp"
#include "MaterialQuad.hpp"
#include "LinPardiso.hpp"
#include "Quadrature.hpp"
#include "StepperNewmark.hpp"
#include "StepperRigid.hpp"
#include "StepperStatic.hpp"
#include "StretchMesh.hpp"

void initLinSolvers(World * world)
{
  std::vector<int> I, J;
  for (size_t i = 0; i < world->em_.size(); i++){
    ElementMesh * em = world->em_[i];
    if (world->solvers.size() <= i){
      LinPardiso<double>  * solver = new LinPardiso<double>();
      solver->init();
      em->stiffnessPattern(I, J);
      Eigen::SparseMatrix<double> K = em->Kpattern;
      solver->init(K);
      world->solvers.push_back(solver);
    }
  }
}

int loadScene(World * world, const ConfigFile & conf)
{
  int dim = 3;
  conf.getInt("dim", dim);
  ElementMesh * em = new ElementMesh();
  em->dim = dim;

  std::string quadratureType = "Gauss2";
  if (conf.hasOpt("quadrature")){
    quadratureType = conf.getString("quadrature");
  }
  const Quadrature* quadrature = 0;
  if (quadratureType == "Gauss2"){
    quadrature = &(Quadrature::Gauss2);
  }
  else if (quadratureType == "NinePt"){
    quadrature = &(Quadrature::NinePt);
  }

  if(dim == 1){
    quadrature = &(Quadrature::Gauss2_2D);
  }

  if (conf.hasOpt("materialfile")){
    std::string materialFile = conf.dir + conf.getString("materialfile");
    FileUtilIn matIn(materialFile);
    if (!matIn.good()){
      return -1;
    }
    std::vector<double> densities;
    std::vector<StrainEne*> ene = loadMaterials(matIn.in, densities);
    matIn.close();
    std::vector<Material *> material(ene.size());
    for (unsigned int ii = 0; ii < material.size(); ii++){
      MaterialQuad * matq = new MaterialQuad();
      if (dim == 2){
        matq->q = quadrature;
        matq->e.resize(matq->q->x.size());
      }
      for (unsigned int jj = 0; jj < matq->e.size(); jj++){
        matq->e[jj] = ene[ii];
      }
      material[ii] = matq;
      material[ii]->density = densities[ii];
    }
    world->materials = material;
  }

  for (size_t i = 0; i < world->materials.size(); i++){
    em->addMaterial(world->materials[i]);
  }

  bool fixbottom = true;
  conf.getBool("fixbottom", fixbottom);
  
  Eigen::Vector3d ff(0, 0, 0);

  if(conf.hasOpt("force")){
    std::vector<float> confForce = conf.getFloatVector("force");
    if(confForce.size() != 3){
      std::cout<<"Wrong force dimension in config file\n";
    }
    for(int ii =0 ; ii<3; ii++){
      ff[ii] = confForce[ii];
    }
  }

  if(!conf.hasOpt("meshfile")){
    std::cout<<"Need a mesh file\n";
    return -1;
  }
  std::string meshfile = conf.dir + conf.getString("meshfile");
  FileUtilIn in(meshfile);
  if(!in.good()){
    return -1;
  }

  float meshscale = 1;
  conf.getFloat("meshscale", meshscale); 
  em->load(in.in,meshscale);

  //apply some forces
  Eigen::Vector3d max(0,0,0), min(100,100,100);
  BBox(em->X, min, max);

  float eps = 1e-5f;
  //fixed vertices
  if(fixbottom){
    for(unsigned int ii =0 ; ii<em->x.size(); ii++){
      if (em->X[ii][1] < min[1] + eps){
        //&& em->X[ii][0] < min[0] + 0.036){
        for(int jj = 0; jj<dim; jj++){
          em->fixedDof[dim * ii + jj] = 1;
        }
      }
    }
  }

  ////only look at collisions of boundary vertices
  bool reduceCollision = true;
  if (dim == 3 && reduceCollision){
    std::vector<int> nnbr(em->x.size(), 0);
    for (unsigned int ii = 0; ii<em->e.size(); ii++){
      Element * ele = em->e[ii];
      for (int jj = 0; jj<ele->nV(); jj++){
        nnbr[ele->at(jj)]++;
      }
    }
    for (unsigned int ii = 0; ii<em->x.size(); ii++){
      if (em->X[ii][1] < min[1] + eps &&
        nnbr[ii]>2//){
        //use full contact for front 1cm vertices.
        && em->X[ii][0]>min[0] + 0.01){
        em->lb[ii][1] = -1e-3;
      }
    }
  }

  //in inches. Distance from the right.
  float forcedist = 0.5;
  float forceHeight = max[1];
  float eleSize = em->X[1][2] - em->X[0][2];
  conf.getFloat("forcedist", forcedist);
  conf.getFloat("forceHeight", forceHeight);
  //can contain duplicates
  std::vector<int> forceVertices;
  std::vector<int> bottomVertices;
  for(unsigned int ii =0 ; ii<em->e.size(); ii++){
    for(int jj =0 ; jj<em->e[ii]->nV(); jj++){
      int vi = em->e[ii]->at(jj);
      Eigen::Vector3d x = em->X[vi];
      if(  forceHeight - x[1]      < eleSize &&
           x[1] - forceHeight < eleSize 
           && (x[0] - min[0]) < forcedist){
        //std::cout << "force v " << vi << "\n";
        forceVertices.push_back(vi);
      }
      if( x[1] - min[1] < eps
          &&(x[0] - min[0])<forcedist){
        bottomVertices.push_back(vi);
      }
    }
  }

  float forceScale = 1.0f/forceVertices.size();
  for(unsigned int ii = 0; ii<forceVertices.size(); ii++){
    int vi = forceVertices[ii];
    em->fe[vi] += forceScale * ff;
  }

  conf.getFloat("dt", world->flight_h);
  conf.getFloat("dtcontact", world->contact_h);
  conf.getBool("handle_contact", world->handle_contact);
  conf.getInt("contactmodel", world->contactModel);
  conf.getFloat("damp_alpha", em->dampAlpha);
  conf.getFloat("damp_beta", em->dampBeta);
  bool specialvert = false;
  conf.getBool("specialvert", specialvert);
  if (specialvert){
    //different lower bound for the special vertex, which is the height of the feet.
    float lb = 0.006f;
    conf.getFloat("speciallb", lb);
    specialVertexBound(em, lb);
  }

  em->check();

  std::vector<Stepper*> stepperList;
  std::vector<std::string> stepperTypes = conf.getStringVector("stepper");

  for(unsigned int ii = 0; ii<stepperTypes.size(); ii++){
    Stepper * stepper = 0;
    if (stepperTypes[ii]=="newmark"){
      StepperNewmark * sn = new StepperNewmark();
      stepper = sn;
      sn->useOldContact = true;
      //sn->useOldContact = false;
      int val=0;
      conf.getInt("frictionModel", val);
      sn->frictionModel = (StepperNewmark::FrictionModel)val;
    }else if(stepperTypes[ii] =="static"){
      StepperStatic * st = new StepperStatic();
      stepper = st;
    }else if(stepperTypes[ii] == "rigid"){
      StepperRigid * st = new StepperRigid();
      stepper = st;
    }
    if(stepper != 0){
      stepperList.push_back(stepper);
      conf.getBool("addgravity", stepper->addGravity);
      stepper->logFileName = conf.getString("logfile");
      conf.getInt("nSteps", stepper->nSteps);
    }
  }

  if(stepperList.size()==0){
    std::cout<<"Error in config: need a time stepper\n";
    return -1;
  }

  if(conf.hasOpt("box")){
    std::vector<float> boxInput = conf.getFloatVector("box");
    for(int dim = 0; dim<3; dim++){
        world->box.mn[dim] = boxInput[dim];
        world->box.mx[dim] = boxInput[3 + dim];
    }
    //world->box.computeEdgeVerts();
  }
  else{
    world->box.mn = Eigen::Vector3d(-10, -10, -10);
    world->box.mx = Eigen::Vector3d(-10, -10, -10);
  }
  
  conf.getBool("terminateOnCollision", world->terminateOnCollision);
  conf.getBool("savepos", world->savePos);
  world->simDirName = conf.dir + "/" + conf.getString("resultdir");
  conf.getInt("frameskip", world->frameSkip);

  if(conf.hasOpt("restpose")){
    std::string filename = conf.getString("restpose");
    std::ifstream in(filename);
    if(!in.good()){
      std::cout<<"Warning: Can't read "<<filename<<"\n";
    }else{
      load_arr(em->X,in);
      em->x = em->X;
    }
    in.close();
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
    if (em->dim == 3){
      segmentJumper3d(*em, opts);
      stretchJumper3d(*em, opts);
    }
  }

  if (conf.hasOpt("initpos")){
    std::string filename = conf.dir + "/" + conf.getString("initpos");
    std::ifstream in(filename);
    if (!in.good()){
      std::cout << "Warning: Can't read " << filename << "\n";
    }
    else{
      load_arr(em->x, in);
      load_arr(em->v, in);
    }
    in.close();
  }
  
  world->quadrature = quadrature;
  em->initElements(quadrature);
  world->em_.push_back(em);

  //for collision detection
  if (world->handle_contact){
    TrigMesh * tm = new TrigMesh();
    hexToTrigMesh(em, tm);
    world->trigm.push_back(tm);
  }

  initLinSolvers(world);
  world->stepperList = stepperList;
  std::cout<<stepperList.size()<<" steppers.\n";

  return 0;
}
