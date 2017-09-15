#include "WorldCompress.hpp"

#include "ArrayUtil.hpp"
#include "ElementMesh.hpp"
#include "ElementMeshUtil.hpp"
#include "FileUtil.hpp"
#include "LinPardiso.hpp"
#include "StepperGrad.hpp"
#include "StepperStatic.hpp"

int loadCompressScene(WorldCompress * w, const ConfigFile & conf)
{
  w->quadrature = &Quadrature::Gauss2;
  loadMaterials(w->materials, conf, w->quadrature);

  if (!conf.hasOpt("meshfile")){
    std::cout << "Need a mesh file\n";
    return -1;
  }
  std::string meshfile = conf.dir + conf.getString("meshfile");
  FileUtilIn in(meshfile);
  if (!in.good()){
    return -1;
  }
  float meshscale = 1;
  conf.getFloat("meshscale", meshscale);
  w->em_ = new ElementMesh();
  w->em_->dim = 3;
  w->em_->load(in.in, meshscale);
  in.close();
  w->em_->initArrays();
  for (size_t j = 0; j < w->materials.size(); j++){
    w->em_->addMaterial(w->materials[j]);
  }

  //fix top and bottom vertices.
  int dim = w->em_->dim;
  Eigen::Vector3d max(0, 0, 0), min(100, 100, 100);
  BBox(w->em_->X, min, max);
  Eigen::Vector3d middle = 0.5 * (max + min);
  float eleSize = (float)(w->em_->X[1][2] - w->em_->X[0][2]);
  float eps = 0.1f * eleSize;
  for (unsigned int ii = 0; ii < w->em_->x.size(); ii++){
    //fix bottom 
    if (w->em_->X[ii][1] < min[1] + eps){
      w->em_->fixedDof[dim * ii + 1] = 1;
    }
    //fix top
    if (w->em_->X[ii][1] > max[1] - eps){
      w->em_->fixedDof[dim * ii + 1] = 1;
    }
    //fix x+
    if (w->em_->X[ii][0] > max[0] - eps){
      w->em_->fixedDof[dim * ii + 0] = 1;
    }
    //fix middle slice along z direction.
    //if ( abs(w->em_->X[ii][2] - middle[2]) < eps){
    //  w->em_->fixedDof[dim * ii + 2] = 1;
    //}
    //fix z+
    if ( w->em_->X[ii][2] > max[2] - eps){
      w->em_->fixedDof[dim * ii + 2] = 1;
    }
  }

  conf.getInt("compressSteps", w->compressSteps);
  conf.getFloat("compress", w->compressRatio);

  w->em_->check();
  w->em_->initElements(w->quadrature);

  StepperStatic * st = new StepperStatic();
  st->handle_contact = false;
  st->nSteps = 10;

  w->stepper = st;
  w->simDirName = conf.dir + "/" + conf.getString("resultdir");
  w->frameCnt = 0;

  bool periodic = false;
  conf.getBool("periodic", periodic);
  if (periodic){
    //add planar constraint for all x=0 and z = 0 vertices.
    std::vector<int> leftDof, backDof;
    for (unsigned int ii = 0; ii < w->em_->x.size(); ii++){
      if (w->em_->X[ii][0] < min[0] + eps){
        leftDof.push_back(dim * ii);
      }
      if (w->em_->X[ii][2] < min[2] + eps){
        backDof.push_back(dim * ii + 2);
      }
    }
    if (leftDof.size() > 1){
      st->CEq.push_back(leftDof);
    }
    if (backDof.size() > 1){
      st->CEq.push_back(backDof);
    }
  }

  w->initSolver();
  st->solver = w->solver;
  conf.getInt("nSteps", st->nSteps);
  conf.getBool("savepos", w->savePos);
  return 0;
}

void WorldCompress::initSolver()
{
  LinPardiso<double> * s = new LinPardiso<double>();
  s->init();
  Eigen::SparseMatrix<double> K = em_->getStiffnessSparse();
  for (int i = 0; i < stepper->CEq.size(); i++){
    equalityConstraint(stepper->CEq[i], em_, K);
  }
  s->init(K);
  solver = s;
}

void WorldCompress::loop()
{
  std::vector<Eigen::Vector3d> x0 = em_->x;
  int dim = em_->dim;
  double totalTime = 0;
  stepper->m = em_;
  stepper->init();
  StepperStatic *staticSolver = (StepperStatic*)(stepper);
  staticSolver->fixed0 = em_->fixedDof;
  int ret = 0;
  
  stepper->solver = solver;
  stepper->m = em_;
  saveMeshState();
  frameCnt++;
  //std::ofstream outf("L.txt");
  for (int cstep = 0; cstep < compressSteps; cstep++){
    std::vector<float> scale(3, 1);
    scale[1] = (1-(cstep+1)*compressRatio)/(1-cstep * compressRatio);
    scaleMeshx(em_, scale);
    stepper->init();
    for (int iter = 0; iter < stepper->nSteps; iter++){
      ret = stepper->stepWrapper();
      if (ret != 0){
        break;
      }
    }
    if (savePos){
      saveMeshState();
    }
    frameCnt++;

  }
  std::cout << "compression status " << ret;
}

void WorldCompress::resetSim()
{
  frameCnt = 0;
  em_->x = em_->X;
  std::fill(em_->v.begin(), em_->v.end(), Eigen::Vector3d::Zero());
  
  //move mesh to starting position
  placeAtOrigin(em_);
  Eigen::Vector3d mn, mx;
  BBox(em_->x, mn, mx);
  ///recompute element Jacobian since rest shape may have changed.
  em_->initElements(quadrature);
}

void WorldCompress::saveMeshState()
{
  std::string prefix = simDirName + "/0m";
  std::string filename = sequenceFilename(frameCnt, 0, prefix.c_str());
  std::ofstream out(filename);
  if (!out.good()){
    //std::cerr<<"cannot open output file "<<filename<<"\n";
    return;
  }
  else{
    saveSimState(em_, out);
  }
  out.close();
}

WorldCompress::WorldCompress()
{
  frameCnt = 0;
  savePos = false;
  compressRatio = 0.005f;
  compressSteps = 10;
}

WorldCompress::~WorldCompress()
{

}