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

Logger * logger = 0;

int main(int argc, char* argv[])
{
  std::string a1;

  if (argc < 2) {
	  return 0;
  }else {
    a1 = argv[1];
  }
  
  if (a1 == "3"){
	//converts an .txt fem mesh to 2 obj surface meshes.
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
	//subdivide an fem mesh uniformly.
	//also saves a surf mesh for preview.
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
	for (size_t i = 0; i < ccd.mdl->vflist.size(); i++) {
		std::cout << ccd.mdl->vflist[i].f << " " << ccd.mdl->vflist[i].v << "\n";
	}
  }
  else if (a1 == "6"){
    //save modal meshes.
	//modal vectors are stored in "modefile".
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
