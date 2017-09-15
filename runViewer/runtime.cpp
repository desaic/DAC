#include "ArrayUtil.hpp"
#include "ConfigFile.hpp"
#include "TrigMesh.hpp"
#include "FileUtil.hpp"
#include "World.hpp"
#include "ElementMesh.hpp"
#include "ElementMeshUtil.hpp"
#include "Render.hpp"
#include "Replay.hpp"
#include "StrainEne.hpp"
#include "StretchMesh.hpp"
#include "Quadrature.hpp"
#include "MaterialQuad.hpp"

#include "ElementRegGrid.hpp"

#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>

int loadReplayScene(Replay * replay, const ConfigFile & conf);

int main(int argc, char* argv[])
{
  const char * filename = "config.txt";
  if(argc>1){
    filename = argv[1];
  }

  ConfigFile conf;
  int status = conf.load(filename);
  if(status<0){
    return -1;
  }

  Replay * replay = new Replay();
  //load replay scene
  status = loadReplayScene(replay, conf);
  if(status<0){
    return -1;
  }
  
  Render render;

  if(conf.hasOpt("screenfile")){
    render.screenFileName = conf.getString("screenfile");
  }
  conf.getInt("mode", render.mode);
  if (render.mode == 1){
    replay->viewStress = true;
  }
  else{
    replay->viewStress = false;
  }
  if (conf.hasOpt("skipframe")){
    conf.getInt("skipframe", replay->skipFrame);
  }
  conf.getBool("capturescreen", replay->captureScreen);
  conf.getBool("loadDisplacement", replay->loadDisplacement);

  render.init(replay);
  replay->launchThread();
  render.setCamera(replay->em[0]);
  render.loop();
  return 0;
}

int loadReplayScene(Replay * replay, const ConfigFile & conf)
{
  int dim = 3;
  conf.getInt("dim", dim);
  
  //make materials
  std::string materialFile = "../../models/material2.txt";
  if (conf.hasOpt("materialfile")){
    materialFile = conf.dir + conf.getString("materialfile");
  }

  replay->repeat = false;
  conf.getBool(conf.getString("repeat"), replay->repeat);
  
  FileUtilIn matIn(materialFile);
  if (!matIn.good()){
    return -1;
  }
  std::vector<double> densities;
  std::vector<StrainEne*> ene = loadMaterials(matIn.in,densities);
  matIn.close();
  const Quadrature* quadrature = &(Quadrature::Gauss2);
  if (dim == 1){
    quadrature = &(Quadrature::Gauss2_2D);
  }
  std::vector<Material *> material(ene.size());
  for (unsigned int ii = 0; ii<material.size(); ii++){
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
 
  if (!conf.hasOpt("meshfile")){
    std::cout << "Need a mesh file\n";
    return -1;
  }

  std::vector<std::string> meshfiles = conf.getStringVector("meshfile");
  //close enough to boundary
  float eps = 1e-5;
  for (size_t i = 0; i < meshfiles.size(); i++){
    std::string meshfile = conf.dir + meshfiles[i];
    FileUtilIn in(meshfile);
    if (!in.good()){
      continue;
    }
    ElementMesh * em = new ElementMesh();
    em->dim = dim;
    em->load(in.in);
    in.close();
    for (size_t j = 0; j < material.size(); j++){
      em->addMaterial(material[j]);
    }
    //apply some forces
    Eigen::Vector3d max(0, 0, 0), min(100, 100, 100);
    BBox(em->X, min, max);
    conf.getFloat("damp_alpha", em->dampAlpha);
    conf.getFloat("damp_beta", em->dampBeta);
    conf.getFloat("force_draw_scale", em->forceDrawingScale);
   
    //temporary fix to ignore coarse fem sim material 
    for (size_t j = 0; j < em->me.size(); j++){
      if (em->me[j] >= material.size()){
        em->me[j] = material.size() - 1;
      }
    }
    em->check();

    replay->em.push_back(em);
    TrigMesh * tm = new TrigMesh();
    hexToTrigMesh(em, tm);
    //if (repeat)

    replay->trigm.push_back(tm);
    //tm->save_obj("mat0.obj", 0, em);
    //tm->save_obj("mat1.obj", 1, em);
    em->initElements(quadrature);
    em->computeMass();
  }

  //in inches. Distance from the right.
  float forcedist = 0.5;
  conf.getFloat("forcedist", forcedist);
  conf.getFloat("dt", replay->flight_h);
  conf.getFloat("dtcontact", replay->contact_h);
  conf.getInt("nSteps", replay->nFrames);
  
  if (conf.hasOpt("box")){
    std::vector<float> boxInput = conf.getFloatVector("box");
    for (int dim = 0; dim<3; dim++){
      replay->box.mn[dim] = boxInput[dim];
      replay->box.mx[dim] = boxInput[3 + dim];
    }
  }

  bool specialvert = false;
  conf.getBool("specialvert", specialvert);
  if (specialvert){
    //different lower bound for the special vertex, which is the height of the feet.
    float lb = 0.006f;
    conf.getFloat("speciallb", lb);
    specialVertexBound(replay->em[0], lb);
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
    if (replay->em[0]->dim == 3){
      segmentJumper3d(*replay->em[0], opts);
      stretchJumper3d(*replay->em[0], opts);
      replay->em[0]->initElements(quadrature);
    }
  }
  //FileUtilOut out("femmesh_opt.txt");
  //replay->em[0]->saveMesh(out.out);
  //out.out.close();
  if (conf.hasOpt("resultdir")){
    replay->simDirName = conf.dir + conf.getString("resultdir");
    //replay->simDirName = conf.getString("resultdir");
  }

  replay->stairHeightStart = -1;
  replay->stairWidth = -1;

  if (conf.hasOpt("stairHeightStart")){
    conf.getFloat("stairHeightStart", replay->stairHeightStart);
    conf.getFloat("stairHeight", replay->stairHeight);
    conf.getFloat("stairWidth", replay->stairWidth);
  }

  replay->quadrature = quadrature;
  conf.getBool("saveObj", replay->saveObj);
  replay->objDir = conf.dir + "/" + conf.getString("objDir");

  conf.getInt("frame0", replay->frame0);
  return 0;
}

