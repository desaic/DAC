/*
 * Render.cpp
 *
 *  Created on: Aug 20, 2013
 *      Author: desaic
 */
#include "ArrayUtil.hpp"
#include "Contact.hpp"
#include "Element.hpp"
#include "ElementMesh.hpp"
#include "FileUtil.hpp"
#include "TrigMesh.hpp"
#include "glheader.hpp"
#include "Render.hpp"
#include "Replay.hpp"
#include "Stepper.hpp"

#include "png_io.h"

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <chrono>

static const int DEFAULT_GL_LINE_WIDTH = 2;
const float epsilon = 0.01f;
Render * render;
//mouse state
bool captureMouse = false;

static void error_callback(int error, const char* description)
{
    fputs(description, stderr);
}

static void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
  if(action == GLFW_PRESS){
    switch(key){
    case GLFW_KEY_W:
      render->camera.keyhold[0]=true;
      break;
    case GLFW_KEY_S:
      render->camera.keyhold[1]=true;
      break;
    case GLFW_KEY_A:
      render->camera.keyhold[2]=true;
      break;
    case GLFW_KEY_D:
      render->camera.keyhold[3]=true;
      break;
    case GLFW_KEY_R:
      render->camera.keyhold[4]=true;
      break;
    case GLFW_KEY_F:
      render->camera.keyhold[5]=true;
      break;
    }
  }else if(action==GLFW_RELEASE){
    switch(key){
    case GLFW_KEY_ESCAPE:
      glfwSetWindowShouldClose(window, GL_TRUE);
      break;
    case GLFW_KEY_W:
      render->camera.keyhold[0]=false;
      break;
    case GLFW_KEY_S:
      render->camera.keyhold[1]=false;
      break;
    case GLFW_KEY_A:
      render->camera.keyhold[2]=false;
      break;
    case GLFW_KEY_D:
      render->camera.keyhold[3]=false;
      break;
    case GLFW_KEY_R:
      render->camera.keyhold[4]=false;
      break;
    case GLFW_KEY_F:
      render->camera.keyhold[5]=false;
      break;

    case GLFW_KEY_O:
      render->getReplay()->resetSim();
      break;

    case GLFW_KEY_P:
      render->getReplay()->pauseSim();
      break;
    case GLFW_KEY_B:
      render->getReplay()->stepBackSim();
      break;
    case GLFW_KEY_L:
      render->getReplay()->saveMeshPly();
      break;
    case GLFW_KEY_V:
      render->getReplay()->saveVObj();
      break;
    case GLFW_KEY_C:
      render->contact_idx++;
      if (render->contact_idx >= render->getReplay()->contact.size()){
        render->contact_idx = 0;
      }
      break;

    case GLFW_KEY_LEFT_BRACKET:
      render->getReplay()->singleStepSim();
      break;
    case GLFW_KEY_RIGHT_BRACKET:
      render->getReplay()->continueSim();
      break;
    }
  }
}

void Render::moveCamera(float dt)
{
  Eigen::Vector3f viewDir = camera.at - render->camera.eye;
  Eigen::Vector3f up = camera.up;
  Eigen::Vector3f right = viewDir.cross(up);
  right[1] = 0;
  viewDir[1] = 0;

  if(camera.keyhold[0]){
    camera.eye += viewDir * dt * camSpeed;
    camera.at  += viewDir * dt * camSpeed;
  }
  if(camera.keyhold[1]){
    camera.eye -= viewDir * dt * camSpeed;
    camera.at  -= viewDir * dt * camSpeed;
  }
  if(camera.keyhold[2]){
    camera.eye -= right * dt * camSpeed;
    camera.at  -= right * dt * camSpeed;
  }
  if(camera.keyhold[3]){
    camera.eye += right * dt * camSpeed;
    camera.at  += right * dt * camSpeed;
  }
  if(camera.keyhold[4]){
    if(camera.eye[1]<10){
      camera.eye[1] += dt * camSpeed;
      camera.at[1]  += dt * camSpeed;
    }
  }
  if(camera.keyhold[5]){
    if(camera.eye[1]>0){
      camera.eye[1] -= dt * camSpeed;
      camera.at[1]  -= dt * camSpeed;
    }
  }
}

double xpos0 , ypos0;

void mouseButtonFun(GLFWwindow *window , int button, int action, int mods)
{
  switch(button){
  case GLFW_MOUSE_BUTTON_LEFT:
    if(action == GLFW_RELEASE){
      captureMouse = !captureMouse;
      if(captureMouse){
        int width, height;
        glfwGetFramebufferSize(window, &width, &height);
        glfwSetCursorPos(window,width/2, height/2);
        xpos0 = width/2;
        ypos0 = height/2;
      }
    }
    break;
  }
}

void mousePosFun(GLFWwindow *window , double xpos, double ypos)
{
  if(captureMouse){
    int width, height;
    glfwGetFramebufferSize(window, &width, &height);
    double dx =  xpos - xpos0;
    double dy = -ypos + ypos0;
    xpos0 = xpos;
    ypos0 = ypos;
//    std::cout<<dx<<" "<<dy<<"\n";
//    glfwSetCursorPos(window,width/2, height/2);
    render->camera.angle_xz += (float)(dx * render->xRotSpeed);
    render->camera.angle_y  += (float)(dy * render->yRotSpeed);
    render->camera.update();
  }
}

void Render::drawEle(int eidx, ElementMesh * m)
{
  Element * ele = m->e[eidx];
  float eleScale = 0.9f;
  Eigen::Vector3d center = Eigen::Vector3d::Zero();
  for(int ii = 0; ii<ele->nV(); ii++){
    center += m->x[ele->at(ii)];
  }
  center/=ele->nV();

  std::vector<std::array<int , 2> > edges = ele->getEdges();
  Eigen::Vector3f col = m->materialColor[m->me[eidx]];
  //col[0] += m->depthColorScale * (float)m->depth[eidx];
  glColor3f(col[0], col[1], col[2]);

  glDisable(GL_LIGHTING);
  glBegin(GL_LINES);
  for(unsigned int ii = 0;ii<edges.size();ii++){
    int vidx = (*ele)[edges[ii][0]];
    if(vidx>=m->x.size()){
      std::cout<<m->x.size()<<"\n";
    }
    Eigen::Vector3d v = m->x[vidx];
    v = eleScale*(v-center)+center;
    glVertex3f((float)v[0], (float)v[1], (float)v[2]);
    v = m->x[(*ele)[edges[ii][1]]];
    v = eleScale*(v-center)+center;
    glVertex3f((float)v[0], (float)v[1], (float)v[2]);
  }

  //if (m->u != 0){
    //glColor3f(col[0], col[1], col[2]);
    //glColor3f(0.2f, 0.7f, 0.1f);
  //  for (unsigned int ii = 0; ii < edges.size(); ii++){
  //    Eigen::Vector3d v = (*(m->u))[eidx][edges[ii][0]];
  //    glVertex3f((float)v[0], (float)v[1], (float)v[2]);
  //    v = (*(m->u))[eidx][edges[ii][1]];
  //    glVertex3f((float)v[0], (float)v[1], (float)v[2]);
  //  }
  //}

  glEnd();
  glEnable(GL_LIGHTING);

}

void Render::drawEleMesh(ElementMesh * eMesh)
{
  for(unsigned int ii = 0;ii<eMesh->e.size();ii++){
    //drawEle(ii,eMesh);
  }
  //draw force
  glColor3f(0.2f, 0.9f, 0.2f);
  glDisable(GL_LIGHTING);
  glBegin(GL_LINES);
  for(unsigned int ii = 0; ii<eMesh->v.size(); ii++){
    Eigen::Vector3d v = eMesh->x[ii];
    glVertex3f((float)v[0], (float)v[1], (float)v[2]);
    v += eMesh->forceDrawingScale * eMesh->fe[ii];
    glVertex3f((float)v[0], (float)v[1], (float)v[2]);
  }

  //visualize contact forces?
  glColor3f(0.5f, 0.5f, 0.2f);
  if(eMesh->fc.size()>0){
    for(unsigned int ii = 0; ii<eMesh->fc.size(); ii++){
      Eigen::Vector3d v = eMesh->x[ii];
      glVertex3f((float)v[0], (float)v[1], (float)v[2]);
      v += eMesh->forceDrawingScale * eMesh->fc[ii];
      //v += eMesh->forceDrawingScale * eMesh->v[ii];
      glVertex3f((float)v[0], (float)v[1], (float)v[2]);
    }
  }

  glEnd();
  glEnable(GL_LIGHTING);
}

void drawTrigMesh(TrigMesh * m, int mode)
{
  if (mode == 0){
    glColor3f(0.3f, 0.4f, 0.9f);
    glDisable(GL_LIGHTING);
    glBegin(GL_LINES);
    const int N_EDGE = 2;
    int edges[N_EDGE][2] = { { 0, 1 }, { 1, 2 }};
    for (unsigned int t = 0; t<m->t.size(); t += 2){

      for (unsigned int ii = 0; ii<N_EDGE; ii++){
        int vidx = m->t[t][edges[ii][0]];
        if (vidx >= m->v.size()){
          std::cout << m->v.size() << "\n";
        }
        Eigen::Vector3f v = m->v[vidx].cast<float>();
        glVertex3f(v[0], v[1], v[2]);
        v = (m->v[m->t[t][edges[ii][1]]]).cast<float>();
        glVertex3f(v[0], v[1], v[2]);
      }
    }
    glEnd();
  }
  else if(mode == 1){
    if (m->tcolor.size() < m->t.size()){
      m->tcolor.resize(m->t.size(), 0.5 * Eigen::Matrix3f::Ones());
    }

    std::vector<Eigen::Vector3f> vcolor(m->v.size(), Eigen::Vector3f::Zero());
    std::vector<int> vcnt(m->v.size(), 0);
    for (size_t i = 0; i < m->t.size(); i++){
      for (int j = 0; j < 3; j++){
        int vidx = m->t[i][j];
        vcolor[vidx] += m->tcolor[i].row(j);
        vcnt[vidx] ++;
      }
    }
    for (size_t i = 0; i < m->v.size(); i++){
      vcolor[i] = (1.0 / vcnt[i]) * vcolor[i];
    }

    GLfloat color[4] = { 0.2f, 0.2f, 0.2f , 1.0f};
    glDisable(GL_LIGHTING);
    glMaterialfv(GL_FRONT, GL_SPECULAR, color);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, color);
    glBegin(GL_TRIANGLES);
    for (unsigned int t = 0; t<m->t.size(); t++){
      Eigen::Vector3f v[3];
      for (unsigned int ii = 0; ii<3; ii++){
        int vidx = m->t[t][ii];
        v[ii] = m->v[vidx].cast<float>();
      }
      Eigen::Vector3f n = (v[1] - v[0]).cross(v[2] - v[0]).normalized();
      
      //color[0] = m->tcolor[t](0, 0);
      //color[1] = m->tcolor[t](0, 1);
      //color[2] = m->tcolor[t](0, 2);
      color[0] = vcolor[m->t[t][0]][0];
      color[1] = vcolor[m->t[t][0]][1];
      color[2] = vcolor[m->t[t][0]][2];
      glNormal3f(n[0], n[1], n[2]);
      glColor3f(color[0], color[1], color[2]);
      glVertex3f(v[0][0], v[0][1], v[0][2]);

      //color[0] = m->tcolor[t](1, 0);
      //color[1] = m->tcolor[t](1, 1);
      //color[2] = m->tcolor[t](1, 2);

      color[0] = vcolor[m->t[t][1]][0];
      color[1] = vcolor[m->t[t][1]][1];
      color[2] = vcolor[m->t[t][1]][2];

      glNormal3f(n[0], n[1], n[2]);
      glColor3f(color[0], color[1], color[2]);
      glVertex3f(v[1][0], v[1][1], v[1][2]);

      //color[0] = m->tcolor[t](2, 0);
      //color[1] = m->tcolor[t](2, 1);
      //color[2] = m->tcolor[t](2, 2);
      color[0] = vcolor[m->t[t][2]][0];
      color[1] = vcolor[m->t[t][2]][1];
      color[2] = vcolor[m->t[t][2]][2];

      glNormal3f(n[0], n[1], n[2]);
      glColor3f(color[0], color[1], color[2]);
      glVertex3f(v[2][0], v[2][1], v[2][2]);
    }
    glEnd();
  }
  else{
    GLfloat color[4] = { 0.6f, 0.7f, 0.9f, 1.0f };
    glMaterialfv(GL_FRONT, GL_SPECULAR, color);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, color);
    glBegin(GL_TRIANGLES);
    for (unsigned int t = 0; t<m->t.size(); t++){
      Eigen::Vector3f v[3];
      for (unsigned int ii = 0; ii<3; ii++){
        int vidx = m->t[t][ii];
#ifdef _DEBUG
        if (vidx >= m->v.size()){
          std::cout << m->v.size() << "\n";
        }
#endif
        v[ii] = m->v[vidx].cast<float>();
      }
      color[2] = t / (float)m->t.size();
      Eigen::Vector3f n = (v[1] - v[0]).cross(v[2] - v[0]).normalized();
      glColor3f(color[0], color[1], color[2]);
      glNormal3f(n[0], n[1], n[2]);
      glVertex3f(v[0][0], v[0][1], v[0][2]);
      glNormal3f(n[0], n[1], n[2]);
      glVertex3f(v[1][0], v[1][1], v[1][2]);
      glNormal3f(n[0], n[1], n[2]);
      glVertex3f(v[2][0], v[2][1], v[2][2]);
    }
    glEnd();
  }
}

void Render::drawContact(Contact * c, const std::vector<ElementMesh * > & em)
{
  glColor3f(0.8f, 0.5f, 0.4f);
  glDisable(GL_LIGHTING);
  if (c->m[1] >= 0){
    glBegin(GL_LINES);
    Eigen::Vector3d v0 = em[c->m[0]]->x[c->v1];
    for (unsigned int i = 0; i<4; i++){
      Eigen::Vector3d v1 = em[c->m[1]]->x[c->v2[i]];
      glVertex3f((GLfloat)v0[0], (GLfloat)v0[1], (GLfloat)v0[2]);
      glVertex3f((GLfloat)v1[0], (GLfloat)v1[1], (GLfloat)v1[2]);
    }
    glEnd();
  }
  
  float len = (float)em[0]->X[1][2] - em[0]->X[0][2];
  int NORMAL_LINE_WIDTH= 8;
  glLineWidth(NORMAL_LINE_WIDTH);
  glBegin(GL_LINES);
  glVertex3f((GLfloat)c->x[0], (GLfloat)c->x[1], (GLfloat)c->x[2]);
  glVertex3f((GLfloat)(c->normal[0] * len + c->x[0]), (GLfloat)(c->normal[1]*len  + c->x[1]), (GLfloat)(c->normal[2] * len + c->x[2]));
  glEnd();
  glLineWidth(DEFAULT_GL_LINE_WIDTH);
  glEnd();

}

void Render::draw()
{
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  gluLookAt(camera.eye[0], camera.eye[1], camera.eye[2],
    camera.at[0], camera.at[1], camera.at[2],
    camera.up[0],camera.up[1],camera.up[2]);

  std::unique_lock<std::mutex> lock(replay->mtx, std::defer_lock);
  lock.lock();

  for (unsigned int ii = 0; ii<replay->em.size(); ii++){
    drawEleMesh(replay->em[ii]);
  }

  for (unsigned int ii = 0; ii<replay->trigm.size(); ii++){
    drawTrigMesh(replay->trigm[ii], mode);
  }

  if (contact_idx < replay->contact.size()){
    drawContact(&(replay->contact[contact_idx]), replay->em);
  }
  lock.unlock();

  //line for axis
  glDisable(GL_LIGHTING);
  glBegin(GL_LINES);
  glColor3f(0,0,0);
  glVertex3f( -10.0f, 0, 0.0f);
  glVertex3f(  10.0f, 0, 0.0f );

  glVertex3f(0.0f, 0, -10.0f);
  glVertex3f(0.0f, 0, 10.0f);

  glVertex3f(0.0f, -10, 0.0f);
  glVertex3f(0.0f, 10, 0.0f);

  //glVertex3f((float)replay->box.mn[0], (float)replay->box.mn[1], 0);
  //glVertex3f((float)replay->box.mn[0], (float)replay->box.mx[1], 0);

  //glVertex3f((float)replay->box.mn[0], (float)replay->box.mx[1], 0);
  //glVertex3f((float)replay->box.mx[0], (float)replay->box.mx[1], 0);

  //glVertex3f((float)replay->box.mx[0], (float)replay->box.mx[1], 0);
  //glVertex3f((float)replay->box.mx[0], (float)replay->box.mn[1], 0);
  //
  //if (replay->stairWidth > 0){
  //  //int nSteps = (int)(replay->stairHeightStart / replay->stairHeight);
  //  //for (int i = 0; i < nSteps; i++){
  //    float y = replay->stairHeightStart;
  //    float x = replay->stairWidth / replay->stairHeight * replay->stairHeightStart;
  //    glVertex3f(0,  y,  0);
  //    glVertex3f(x, 0, 0);
  //  //}
  //}
  //
  
  glEnd();

  GLfloat floorCol[4] = { 1, 1, 1, 1 };
  glEnable(GL_LIGHTING);

  GLfloat position[] = { 2.0f, 2, 2, 1.0f };
  GLfloat position1[] = { -1.0f, -1, -1, 1.0f };
  glLightfv(GL_LIGHT0, GL_POSITION, position);
  glLightfv(GL_LIGHT1, GL_POSITION, position1);

  glMaterialfv(GL_FRONT, GL_SPECULAR, floorCol);
  glMaterialfv(GL_FRONT, GL_DIFFUSE, floorCol);
  GLfloat s = 10;
  glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, &s);

  //glBegin(GL_TRIANGLE_STRIP);
  //glNormal3f(0, 1, 0);
  //glVertex3f( -1.0f, -.5f , -1.0f);
  //glVertex3f( -1.0f, -.5f,  1.0f );
  //glVertex3f(  1.0f, -.5f , -1.0f);
  //glVertex3f(  1.0f, -.5f ,  1.0f);
  //glEnd();

  glfwSwapBuffers(window);
}

int Render::loop()
{
  std::chrono::milliseconds::rep t0, t1, dt;
  auto duration = std::chrono::system_clock::now().time_since_epoch();
  t0 = std::chrono::duration_cast<std::chrono::milliseconds>(duration).count();

  while (!glfwWindowShouldClose(window))
  {
    draw();
    glfwPollEvents();

    duration = std::chrono::system_clock::now().time_since_epoch();
    t1 = std::chrono::duration_cast<std::chrono::milliseconds>(duration).count();
    dt = t1 - t0;

    moveCamera((float)dt);
    t0 = t1;
    if(replay->captureScreen){
      screenshot();
      replay->captureScreen = false;
    }
    if(replay->state == Replay::DONE){
      if(replay->em.size()>0){
        setCamera(replay->em[0]);
        draw();
        draw();
      }
      break;
    }

  }
  glfwDestroyWindow(window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}

void Render::screenshot()
{
  int height = 0;
  int width = 0;

  glfwGetWindowSize(window, &width, &height);

  GLubyte* pixels = new GLubyte[ 3 * width * height];
  int colorFlag=0;
#ifdef __APPLE__
  colorFlag = GL_BGR;
#else
  colorFlag = GL_BGR_EXT;
#endif
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
  glPixelStorei(GL_PACK_ALIGNMENT, 1);
  std::cout << "img size " << width << " " << height << "\n";
  glReadPixels(0, 0, width, height, colorFlag, GL_UNSIGNED_BYTE, pixels);
  screenFileName = sequenceFilename(screenIdx, 0, "screen", 4, ".png");
  save_png(screenFileName, pixels, width, height);

  screenIdx++;
}

void window_size(GLFWwindow* window, int width, int height)
{
  float ratio;
  ratio = width / (float) height;
  glViewport(0, 0, width, height);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(40.0, ratio, 0.01, 20.0);
}

void
Render::init(Replay * _replay)
{
  glfwSetErrorCallback(error_callback);
  if (!glfwInit()){
    std::cout<<"Cannot init glfw\n";
    return;
  }
  glfwWindowHint(GLFW_DEPTH_BITS,32);
  window = glfwCreateWindow(800, 600, "Fem Viewer", NULL, NULL);
  if (!window)
  {
    glfwTerminate();
    std::cout<<"Cannot create window\n";
    return;
  }
  replay = _replay;
  render = this;
  glfwMakeContextCurrent(window);
  glfwSetKeyCallback(window, key_callback);
  glfwSetCursorPosCallback( window,mousePosFun );
  glfwSetMouseButtonCallback( window,mouseButtonFun );
  glfwSetWindowSizeCallback(window, window_size);

  float ratio;
  int width, height;
  glfwGetFramebufferSize(window, &width, &height);
  ratio = width / (float) height;
  glViewport(0, 0, width, height);
  glClearColor(1.0f, 1.0f, 1.0f, 0.0f);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(40.0, ratio, 0.01, 20.0);

  glLineWidth(DEFAULT_GL_LINE_WIDTH);

  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_LIGHT1);
  glEnable(GL_DEPTH_TEST);

  GLfloat white[]={1.0f, 1.0f, 1.0f, 1.0f};
  GLfloat grey[] ={0.3f, 0.3f, 0.3f, 1.0f};

  glLightfv (GL_LIGHT1, GL_DIFFUSE, white);
  glLightfv (GL_LIGHT1, GL_SPECULAR, white);
  glLightfv (GL_LIGHT0, GL_AMBIENT, grey);
}

void
Render::setCamera(ElementMesh * m)
{
  Eigen::Vector3d mn, mx;
  BBox(m->x, mn, mx);
  Eigen::Vector3d center, size;
  center=0.5f*(mx + mn);
  size = mx - mn;
  camera.at = center.cast<float>();
  float maxsize = (float)std::max(size[0], size[1]);
  camera.eye = camera.at;
  camera.eye[2] += 1.5f*maxsize + 0.5f*(float)size[2];
  camera.update();
  camSpeed = std::abs(camera.eye[2]) * 5e-4;
}

void
Render::setCamera(TrigMesh * m)
{
  Vector3s mn, mx;
  BBox( m->v, mn, mx);
  Eigen::Vector3f center, size;
  center=0.5*(mx + mn).cast<float>();
  size = (mx - mn).cast<float>();
  camera.at = center.cast<float>();
  float maxsize = std::max(size[0], size[1]);
  camera.eye = camera.at;
  camera.eye[2] += 1.5f*maxsize + 0.5f*size[2];
  camera.update();
  camSpeed = std::abs(camera.eye[2]) * 5e-4;
}

Render::Render():replay(0),anim(false),mode(1),
  xRotSpeed(0.004f),
  yRotSpeed(0.004f),camSpeed(0.00005f),
  contact_idx(0),
  screenIdx(0),
  screenFileName ("screen.png")
{}

Render::~Render()
{}
