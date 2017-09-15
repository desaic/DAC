#ifndef REPLAY_HPP
#define REPLAY_HPP
#include <thread>
#include <mutex>
#include <condition_variable>
#include <vector>

#include "Contact.hpp"
#include "Quadrature.hpp"
#include "World.hpp"
#include "ElementMesh.hpp"
#include "TrigMesh.hpp"

class Replay
{
public:
  enum State{ PAUSE, SINGLE, ALL, DONE };
  State state;
  
  std::vector<ElementMesh * > em;
  std::vector<TrigMesh* > trigm;
  std::vector<Contact> contact;
  Box box;
  std::thread th;
  std::mutex mtx;
  std::condition_variable cv;

  std::string   simDirName;
  int padding;
  float flight_h, contact_h;
  const Quadrature * quadrature;
  int frame0;
  int frame;
  int nFrames;

  bool repeat;

  ///@brief whether save screenshots
  bool captureScreen;
  bool viewStress;
  bool saveObj;

  int skipFrame;
  std::string objDir;
  ///@brief loads displacement instead of positions.
  bool loadDisplacement;
  Replay() :state(PAUSE), padding(0),
    flight_h(1e-3f), contact_h(1e-3f), frame0(0),frame(0),
    nFrames(0), captureScreen(false), viewStress(false), saveObj(false),
    skipFrame(1),
    loadDisplacement(false), stairWidth(-1)
  {}

  ///@param startIdx start looking for mesh at index.
  ///@param meshIdx which mesh in em to update.
  int loadMesh(int startIdx, int meshIdx);
  void run();
  void launchThread();
  void pauseSim();
  void singleStepSim();
  void continueSim();
  void stepBackSim();
  void resetSim();
  void saveMeshPly();
  void saveVObj();
  float stairWidth;
  float stairHeight;
  float stairHeightStart;

};

#endif
