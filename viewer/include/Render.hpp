/*
 * Render.hpp
 *
 *  Created on: Aug 20, 2013
 *      Author: desaic
 */

#ifndef RENDER_HPP_
#define RENDER_HPP_

#include <vector>
#include <Eigen/Dense>

#include "Camera.hpp"
struct GLFWwindow;
class ElementHex;
class Replay;
class Element;
class ElementMesh;
class TrigMesh;
class Stepper;
struct Contact;

class Render
{
public:
  Render();
  void init(Replay * replay);
  int loop();
  //set camera to look at the mesh
  void setCamera(ElementMesh * m);
  void setCamera(TrigMesh * m);
  void draw();
  void drawEle(int eidx, ElementMesh * eMesh);
  void drawEleMesh(ElementMesh * eMesh);
  void drawContact(Contact * c, const std::vector<ElementMesh * > & em);

  void moveCamera(float dt);
  Replay * getReplay(){ return replay; }
  virtual ~Render();

  void screenshot();

  bool anim;
  //0 triangle mesh mode
  //1 default. Wireframe mode.
  //2 glColor. display stress.
  int mode;
  Camera camera;
  GLFWwindow* window;
  
  std::vector<Eigen::Vector3f> matColor;

  //how fast to rotate in x and y axis
  float xRotSpeed, yRotSpeed;
  float camSpeed;
  //which contact to draw.
  int contact_idx;
  int screenIdx;
  std::string screenFileName;
private:
  ///@brief Render does not own this pointer.
  Replay * replay;
};

#endif /* RENDER_HPP_ */
