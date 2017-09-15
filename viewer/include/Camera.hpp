#ifndef CAMERA_HPP
#define CAMERA_HPP

#include <math.h>
#include <Eigen/Dense>

const float MAXY = 0.49f * 3.141593f;
const float MAXXZ = 2*3.141593f;

#define NKEYS 6

struct Camera{
  Camera():angle_xz(3.14f),angle_y(0){
    for (int ii=0;ii<3;ii++){
      eye[ii]=0.0f;
      at[ii]=0.0f;
      up[ii] = 0.0f;
    }
    at[2]=-1.0f;
    at[1] = 0.5f;
    eye[2]= -2.0f;
    eye[1]=0.5f;
    up[1] = 1.0f;

    for(int ii= 0;ii<NKEYS;ii++){
      keyhold[ii] = false;
    }
  }
  float angle_xz, angle_y;
  void update(){
    Eigen::Vector3f viewDir ( -sin(angle_xz), 0, cos(angle_xz));
    if(angle_y > MAXY){
      angle_y = MAXY;
    }
    if(angle_y< -MAXY){
      angle_y = -MAXY;
    }

    if(angle_xz >  MAXXZ){
      angle_xz -=  MAXXZ;
    }
    if(angle_xz < -MAXXZ){
      angle_xz +=  MAXXZ;
    }

    Eigen::Vector3f yaxis = Eigen::Vector3f(0,1,0);
    Eigen::Vector3f right = viewDir.cross(yaxis);
    viewDir = cos(angle_y)*viewDir + sin(angle_y)*yaxis;
    at = eye + viewDir;
    up = right.cross(viewDir);
  }

  Eigen::Vector3f eye;
  Eigen::Vector3f at;
  Eigen::Vector3f up;

  //if keys are held
  //W S A D R F
  bool keyhold[NKEYS];
};

#endif
