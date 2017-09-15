#include "ArrayUtil.hpp"
#include "CCDTrigMesh.hpp"
#include "ConfigFile.hpp"
#include "ConfigUtil.hpp"
#include "ElementHex.hpp"
#include "SampleForces.hpp"
#include "World.hpp"
#include "WorldC.hpp"
#include "WorldCompress.hpp"
#include "WorldStair.hpp"
#include "EigenUtil.hpp"
#include "ElementRegGrid.hpp"
#include "ElementMeshUtil.hpp"
#include "FileUtil.hpp"
#include "logger.h"
#include "MaterialQuad.hpp"
#include "StretchMesh.hpp"
#include "VoxelIO.hpp"
#include "WorldC.hpp"
#include "DeformModel.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>

extern void makeSphere();
void assignBarMat(ElementMesh * em, Eigen::Vector3d center);

Logger * logger = 0;

void genFEMMesh()
{
  //ElementRegGrid grid(10, 10, 1);
  //scaleMesh(&grid, 0.1);
  //std::vector<Element*>e;
  //for (int i = 0; i < 10; i++){
  //  e.push_back(grid.e[i]);
  //}
  //for (int i = 1; i < 10; i++){
  //  e.push_back(grid.e[10 * i]);
  //}
  //for (int i = 1; i < 10; i++){
  //  e.push_back(grid.e[9 + 10 * i]);
  //}
  //grid.e = e;
  //grid.rmEmptyVert();
  //std::ofstream out("L10.txt");
  //grid.saveMesh(out);
  //out.close();

  //ElementRegGrid grid(1, 1, 1);
  //scaleMesh(&grid, 0.1);
  //std::ofstream out("cube.txt");
  //grid.saveMesh(out);
  //out.close();

  //ElementRegGrid grid(1, 6, 1);
  //scaleMesh(&grid, 0.1);
  //std::vector<Element*>e;
  //e.push_back(grid.e[0]);
  //e.push_back(grid.e[2]);
  ////e.push_back(grid.e[4]);
  //grid.e = e;
  //grid.rmEmptyVert();
  //std::ofstream out("2cube.txt");
  //grid.saveMesh(out);
  //out.close();

  //ElementRegGrid grid(2, 3, 1);
  //scaleMesh(&grid, 0.1);
  //std::vector<Element*>e;
  //e.push_back(grid.e[0]);
  //e.push_back(grid.e[1]);
  //e.push_back(grid.e[2]);
  //e.push_back(grid.e[3]);
  //e.push_back(grid.e[5]);

  //grid.e = e;
  //grid.rmEmptyVert();
  //std::ofstream out("c2x3.txt");
  //grid.saveMesh(out);
  //out.close();

  ElementRegGrid grid(4, 4, 4);
  scaleMesh(&grid, 0.03);
  std::ofstream out("cube4.txt");
  grid.saveMesh(out);
  out.close();
}

int main(int argc, char* argv[])
{
  std::string a1;

  if (argc > 1){
    a1 = argv[1];
  }
  //genFEMMesh();
  //return 0;

  //convert .bin microstructures to fem.txt for simulation.
  if (a1 == "1"){
    std::string inname = "tmp.bin";
    std::string outname = "fem.txt";
    if (argc > 2){
      inname = argv[2];
    }
    if (argc > 3){
      outname = argv[3];
    }
    voxel2em(inname, outname);
    return 0;
  }
  else if (a1 == "2"){  
    makeSphere();
    return 0;
  }
  else if (a1 == "3"){
    ElementMesh * em = new ElementMesh();
    std::string inname = "in.txt";
    if (argc > 2){
      inname = argv[2];
    }
    FileUtilIn in(inname);
    em->load(in.in);
    in.close();
    
    savePartObj(em, 0, "surf0.obj");
    savePartObj(em, 1, "surf1.obj");
    return 0;
  }
  else if (a1 == "4"){
    ElementMesh * em = new ElementMesh();
    std::string inname = "in.txt";
    if (argc > 2){
      inname = argv[2];
    }
    FileUtilIn in(inname);
    em->load(in.in);
    in.close();
    ElementMesh * fine=subdivideHex(em, 6);
    fine->initArrays();
    FileUtilOut out("fem6.txt");
    fine->saveMesh(out.out);
    TrigMesh * tm = new TrigMesh;
    hexToTrigMesh(fine, tm);
    tm->save_obj("surf6.obj");
    return 0;
  }
  else if (a1 == "5"){
    //test mccd
    CCDTrigMesh ccd;
    TrigMesh tm[2];
    logger = new Logger("ccdlog.txt");
      for (int i = 0; i < 2; i++) {
        std::string filename(argv[i + 2]);
        std::ifstream in(filename);
        if (!in.good()) {
          std::cout << "Can't open " << filename << "\n";
          return -1;
        }
        tm[i].load(in);
        in.close();
      }
      ccd.init((double *)tm[0].v.data(), tm[0].v.size(), (int*)tm[0].t.data(), tm[0].t.size());
      ccd.update((double *)tm[1].v.data(), tm[1].v.size());
      for (size_t i = 0; i < ccd.mdl->vflist.size(); i++){
        std::cout << ccd.mdl->vflist[i].f << " " << ccd.mdl->vflist[i].v << "\n";
      }
  }
  else if (a1 == "6"){
    //save mode meshes.
    std::string conffile = argv[2];
    ConfigFile conf;
    int status = conf.load(conffile);
    if (status<0){
      return -1;
    }
    if (!conf.hasOpt("meshfile")){
      std::cout << "Need a mesh file\n";
      return -1;
    }
    std::string meshfile = conf.dir + conf.getString("meshfile");
    FileUtilIn in(meshfile);
    if (!in.good()){
      return -1;
    }

    ElementMesh * em = new ElementMesh();
    em->load(in.in, 1.0f);

    std::string modefile = conf.dir + conf.getString("modefile");
    saveModalMeshes(em, modefile);
    return 0;
  }
  else if (a1 == "7"){
    //force-based sampling of a microstructure with non-linear materials.
    std::string filename = "config.txt";
    filename = argv[2];
    ConfigFile conf;
    int status = conf.load(filename);
    if (status<0){
      return -1;
    }
    SampleForces(conf);
    return 0;
  }
  else if (a1 == "8"){
    std::string filename = "config.txt";
    filename = argv[2];
    ConfigFile conf;
    int status = conf.load(filename);
    if (status<0){
      return -1;
    }
    fitParams(conf);
    return 0;
  }
  else if (a1 == "9"){
    std::cout << "extract fine material combinations. "
      <<"Assuming a rectangular grid with "
      <<"dimensions specified in config.txt\n";
    std::string filename = "config.txt";
    filename = argv[2];
    ConfigFile conf;
    int status = conf.load(filename);
    if (status<0){
      return -1;
    }
    extractComb(conf);
    return 0;
  }
  else if (a1 == "10"){
    std::cout << "Simulate coarsened mesh\n";
    std::string filename = "config.txt";
    filename = argv[2];
    ConfigFile conf;
    int status = conf.load(filename);
    if (status<0){
      return -1;
    }
    World * world = new World();
    std::vector<std::vector<double> > params;
    std::string coarseMatFile = conf.dir + "/" + conf.getString("matParam");
    params = loadArr2d<double>(coarseMatFile);
    makeCoarseMats(params, world->materials);
    status = loadScene(world, conf);
    if (status<0){
      return -1;
    }
    //saveStiffness(world->em_[0]);
    Eigen::Vector3d center(0.008, 0.024, 0.008);
    assignBarMat(world->em_[0], center);
    world->loop();
    return 0;
  }

  std::string filename = "config.txt";
  filename = argv[1];
  ConfigFile conf;
  int status = conf.load(filename);
  if(status<0){
    return -1;
  }

  std::string scene = conf.getString("scene");
  if (scene == "jumper"){
    World * world = new World();
    status = loadScene(world, conf);
    if(status<0){
      return -1;
    }

    //world->stage = World::SIM_FLIGHT;
    world->loop();
  }
  else if (scene == "stair"){
    WorldStair * world = new WorldStair();
    status = loadStairScene(world, conf);
    //return 0;
    world->resetSim();
    world->loop();
  }
  else if (scene == "compress"){
    WorldCompress * w = new WorldCompress();
    status = loadCompressScene(w, conf);
    w->loop();
  }
  else{
    WorldC * world = new WorldC();
    loadCScene(world, conf);
    world->resetSim();
    world->loop();
  }

  return 0;
}

void assignBarMat(ElementMesh * em, Eigen::Vector3d center)
{
  for (size_t i = 0; i < em->e.size(); i++){
    Element * ele = em->e[i];
    Eigen::Vector3d eCenter = Eigen::Vector3d::Zero();
    for (int j = 0; j < ele->nV(); j++){
      eCenter += em->X[ele->at(j)];
    }
    eCenter = (1.0 / ele->nV()) * eCenter;
    Eigen::Vector3d disp = eCenter - center;
    int pred1 = (disp[0] < -disp[2]);
    int pred2 = (disp[0] > disp[2]);
    int pred = pred1 + pred2 * 2;
    switch (pred){
    case 0:
      em->me[i] = 1;
      break;
    case 1:
      em->me[i] = 3;
      break;
    case 2:
      em->me[i] = 2;
      break;
    case 3:
    default:
      em->me[i] = 0;
      break;
    }
  }
}