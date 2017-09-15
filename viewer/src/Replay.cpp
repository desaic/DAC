#include "ArrayUtil.hpp"
#include "Contact.hpp"
#include "Element.hpp"
#include "ElementMeshUtil.hpp"
#include "FileUtil.hpp"
#include "Replay.hpp"
#include <string> 
#include <iomanip>

extern const int sw[8][3] =
{ { -1, -1, -1 },
{ -1, -1, 1 },
{ -1, 1, -1 },
{ -1, 1, 1 },
{ 1, -1, -1 },
{ 1, -1, 1 },
{ 1, 1, -1 },
{ 1, 1, 1 }
};

void runReplay(Replay * replay)
{
  replay->run();
}

void Replay::launchThread()
{
  th = std::thread(runReplay, this);
}

void deformEmbed(ElementMesh * em, std::vector<Eigen::Vector3d> & v,
  const std::vector<int> & eidx, const std::vector<Eigen::Vector3d> & natCoord)
{
  for (size_t i = 0; i < v.size(); i++){
    Element * ele = em->e[eidx[i]];
    Eigen::Vector3d n = natCoord[i];
    Eigen::Vector3d newv = Eigen::Vector3d::Zero();
    for (int vi = 0; vi < ele->nV(); vi++){
      float weight = (1.0f / 8) * (1 + sw[vi][0] * n[0])
        *(1 + sw[vi][1] * n[1]) *(1 + sw[vi][2] * n[2]);
      newv += weight * em->x[ele->at(vi)];
    }
    v[i] = newv;
  }
}

void repeatTrigMesh(TrigMesh * tm, Eigen::Vector3d disp, bool mirror_x, bool mirror_y)
{
  //for (int i = 0; i < (int)tm->v.size(); i++){
  //  Eigen::Vector3d newv = tm->v[i];
  //  if (mirror_x)
  //}

  //for (int i = 0; i < (int)tm->t.size(); i++){

  //}
}

int leftTop(const std::vector<Eigen::Vector3d> & v,
  Eigen::Vector3d mn, Eigen::Vector3d mx)
{
  int outIdx = 0;
  //vertices on the left
  std::vector<int> leftIdx;
  float eps = 1e-3;
  for (size_t i = 0; i < v.size(); i++){
    if (v[i][0] < mn[0] + eps){
      leftIdx.push_back(i);
    }
  }
  float maxy = mn[1];
  for (size_t i = 0; i < leftIdx.size(); i++){
    int vidx = leftIdx[i];
    if (v[vidx][1] > maxy){
      maxy = v[vidx][1];
      outIdx = vidx;
    }
  }
  return outIdx;
}

int backTop(const std::vector<Eigen::Vector3d> & v,
  Eigen::Vector3d mn, Eigen::Vector3d mx)
{
  int outIdx = 0;
  //vertices on the left
  std::vector<int> backIdx;
  float eps = 1e-3;
  for (size_t i = 0; i < v.size(); i++){
    if (v[i][2] < mn[2] + eps){
      backIdx.push_back(i);

    }
  }
  float maxy = mn[1];
  for (size_t i = 0; i < backIdx.size(); i++){
    int vidx = backIdx[i];
    if (v[vidx][1] > maxy){
      maxy = v[vidx][1];
      outIdx = vidx;
    }
  }
  return outIdx;
}

int rightBottom(const std::vector<Eigen::Vector3d> & v,
  Eigen::Vector3d mn, Eigen::Vector3d mx)
{
  int outIdx = 0;
  //vertices on the left
  std::vector<int> rightIdx;
  float eps = 1e-3;
  for (size_t i = 0; i < v.size(); i++){
    if (v[i][0] > mx[0] - eps){
      rightIdx.push_back(i);
    }
  }
  float miny = 2 * mx[1];
  for (size_t i = 0; i < rightIdx.size(); i++){
    int vidx = rightIdx[i];
    if (v[vidx][1] < miny){
      miny = v[vidx][1];
      outIdx = vidx;
    }
  }
  return outIdx;
}

int topOrigin(const std::vector<Eigen::Vector3d> & v,
  Eigen::Vector3d mn, Eigen::Vector3d mx)
{
  int outIdx = 0;
  //vertices on the left
  std::vector<int> backIdx;
  float eps = 1e-3;
  for (size_t i = 0; i < v.size(); i++){
    if (v[i][1] > mx[1] - eps){
      backIdx.push_back(i);
    }
  }
  float minDist = 10;;
  for (size_t i = 0; i < backIdx.size(); i++){
    int vidx = backIdx[i];
    Eigen::Vector2d vi;
    vi[0] = v[vidx][0];
    vi[1] = v[vidx][2];
    double dist = vi.squaredNorm();
    if (dist < minDist){
      minDist = dist;
      outIdx = vidx;
    }
  }
  return outIdx;
}

int topRight(const std::vector<Eigen::Vector3d> & v,
  Eigen::Vector3d mn, Eigen::Vector3d mx)
{
  int outIdx = 0;
  //vertices on the left
  std::vector<int> topIdx;
  float eps = 1e-3;
  for (size_t i = 0; i < v.size(); i++){
    if (v[i][1] > mx[1] - 0.01){
      topIdx.push_back(i);
    }
  }
  double mxx = 0;
  for (size_t i = 0; i < topIdx.size(); i++){
    int vidx = topIdx[i];
    double vi = v[vidx][0];
    if (vi > mxx){
      mxx = vi;
      outIdx = vidx;
    }
  }
  return outIdx;
}

int closest(const std::vector<Eigen::Vector3d> & X,
  Eigen::Vector3d target)
{
  double minDist = 1e10;
  int minIdx = 0;
  for (size_t i = 0; i<X.size(); i++){
    double dist = (X[i] - target).squaredNorm();
    if (dist < minDist){
      minDist = dist;
      minIdx = i;
    }
  }
  return minIdx;
}

void Replay::run()
{
  frame = frame0;
  int cnt = 0;
  int objCnt = 0;
  Eigen::Vector3d mn, mx;

  ElementMesh* em0 = em[0];
  
  int e0 = em0->e.size();
  
  BBox(em[0]->X, mn, mx);
  if (saveObj){
    for (int mi = 0; mi < (int)em.size(); mi++){
      for (int mat = 0; mat < 2; mat++){
        std::stringstream objName;
        objName << objDir << "/" << std::to_string(mi) << "m" << mat << "_" << std::setw(4) << std::setfill('0') << std::to_string(objCnt) << ".obj";
        savePartObj(em[mi], mat, objName.str());
      }
    }
    objCnt++;
  }
  while (1){
    std::unique_lock<std::mutex> lck(mtx);
    while (state == PAUSE){
      cv.wait(lck);
    }
    if (state == SINGLE){
      state = PAUSE;
    }
    lck.unlock();
    int status = 0;
    bool loadMore = false;
    for (int mi = 0; mi < (int)em.size(); mi++){
      std::unique_lock<std::mutex> lock(mtx);
      status = loadMesh(frame, mi);
      Eigen::Vector3d p = em[mi]->linearMomentum();
      Eigen::Vector3d center = em[mi]->centerOfMass();
      Eigen::Vector3d L = em[mi]->angularMomentum(center);
      Eigen::Vector3d v = em[mi]->linearVelocity();
      double E = em[mi]->getAllEnergy();
      double ang = em[mi]->getTopAng();
      std::cout << "E " << E << "\n";
      std::cout << "p " << p[0] <<" " <<p[1] <<" " <<p[2] << "\n";
      std::cout << "L " << L[0] << " " << L[1] << " " << L[2] << "\n";
      std::cout << "v " << v[0] << " " << v[1] << " " << v[2] << "\n";
      std::cout << "Height " << center[1] << "\n";

	  if (status >= 0){
        loadMore = true;
        frame = status;
        if (viewStress){
          std::vector<std::vector<double> > P;
          hexToTrigStress(em[mi], trigm[mi], P);
          float maxStress = 0;
          int maxE = 0;
          for (int i = 0; i < P.size(); i++){
            for (int j = 0; j < P[i].size(); j++){
              if (P[i][j] > maxStress){
                maxStress = P[i][j];
                maxE = i;
              }
            }
          }
          std::cout << maxE << " " << maxStress << "\n";
        }
        
        updateSurfVert(em[mi]->x, trigm[mi], trigm[mi]->vidx);
      }
      lock.unlock();
      BBox(em[mi]->x, mn, mx);
      std::cout << "len " << (mx[1] - mn[1]) << "\n";
    }
    if (! loadMore){
      pauseSim();
    }
    
    if (saveObj){ 
      for (int mi = 0; mi < (int)em.size(); mi++){
        for (int mat = 0; mat < 2; mat++){
          std::stringstream objName;
          objName << objDir << "/" << std::to_string(mi) << "m" <<mat<<"_"<< std::setw(4) << std::setfill('0') << std::to_string(objCnt) << ".obj";
          savePartObj(em[mi], mat, objName.str());
        }
      }
      objCnt++;
    }
    cnt++;
    std::string contact_filename = sequenceFilename(frame, simDirName.c_str(),
      "c", padding);
    std::ifstream in(contact_filename);
    std::unique_lock<std::mutex> lock(mtx);
    if (in.good()){
      loadContact(in, contact);
    }
    else{
      contact.clear();
    }
    lock.unlock();
    in.close();
    if (frame < nFrames - 1){
      frame+=skipFrame;
    }
  }

}

void Replay::stepBackSim()
{
  frame-=2;
  if (frame < frame0){
    frame = frame0;
  }
  singleStepSim();

}

void Replay::resetSim()
{
  frame = frame0;
  singleStepSim();
}

void Replay::pauseSim()
{
  std::unique_lock<std::mutex> lck(mtx, std::defer_lock);
  lck.lock();
  state = PAUSE;
  lck.unlock();
}

void Replay::singleStepSim()
{
  std::unique_lock<std::mutex> lck(mtx, std::defer_lock);
  lck.lock();
  if (state == PAUSE){
    state = SINGLE;
    cv.notify_one();
  }
  else{
    state = SINGLE;
    lck.unlock();
  }
}

void Replay::continueSim()
{
  std::unique_lock<std::mutex> lck(mtx, std::defer_lock);
  lck.lock();
  if (state == PAUSE){
    state = ALL;
    cv.notify_one();
  }
  else{
    state = ALL;
    lck.unlock();
  }
}

void Replay::saveMeshPly()
{
  std::string filename = sequenceFilename(frame,0,0,0,".ply");
  if (trigm.size() == 0){
    return;
  }
  savePly(trigm[0], filename);
}

void Replay::saveVObj()
{
  std::string filename = sequenceFilename(frame, 0, 0, 0, ".obj");
  filename = "v" + filename;
  if (em.size() == 0){
    return;
  }
  saveVelObj(em[0], filename);
}

int Replay::loadMesh(int startIdx, int meshIdx)
{
  std::ifstream in;
  int filenum = startIdx;
  std::string filePrefix = std::to_string(meshIdx) + "m";
  nFrames = startIdx + 100;
  while (filenum<nFrames){
    std::string filename = sequenceFilename(filenum, simDirName.c_str(),
      filePrefix.c_str(), padding);
    in.open(filename);
    if (!in.good()){
      filename = sequenceFilename(filenum, simDirName.c_str(),
        filePrefix.c_str(), 0);
      in.open(filename);
    }
    if (!in.good()){
      in.close();
      filenum++;
    }
    else{
      loadSimState(em[meshIdx], in, loadDisplacement);
      in.close();
      break;
    }
  }
  std::cout << filenum << "\n";
  if (filenum >= nFrames){
    //no more files
    return -1;
  }
  return filenum;
}
