#include "Contact.hpp"
#include "Element.hpp"
#include "ElementMesh.hpp"
#include "ElementMeshUtil.hpp"
#include "TrigMesh.hpp"
#include "PtTrigIntersect.hpp"

extern int HexFaces[6][4];

void saveContact(std::ostream & out, const std::vector<Contact> & contact)
{
  out << contact.size() << "\n";
  for (size_t i = 0; i < contact.size(); i++){
    out << contact[i].m[0] << " " << contact[i].m[1] << "\n";
    out << contact[i].v1 << "\n";
    out << contact[i].v2[0] << " " << contact[i].v2[1] << " " << contact[i].v2[2] << " " << contact[i].v2[3] << "\n";
    out << contact[i].x[0] << " " << contact[i].x[1] << " " << contact[i].x[2] << "\n";;
    out << contact[i].normal[0] << " " << contact[i].normal[1] << " " << contact[i].normal[2] << "\n";
  }
}

void loadContact(std::istream & in, std::vector<Contact> & contact)
{
  int N = 0;
  in >> N;
  contact.resize(N);
  for (size_t i = 0; i < contact.size(); i++){
    in >> contact[i].m[0] >> contact[i].m[1];
    in >> contact[i].v1;
    in >> contact[i].v2[0] >> contact[i].v2[1] >> contact[i].v2[2] >> contact[i].v2[3];
    in >> contact[i].x[0] >> contact[i].x[1] >> contact[i].x[2];
    in >> contact[i].normal[0] >> contact[i].normal[1] >> contact[i].normal[2];
  }
}

bool addContact(std::vector<ElementMesh *> & em_,
  std::vector<TrigMesh *> trigm,
  int m0, int m1, int t0, int v1,
  const std::vector<std::vector<Eigen::Vector3d> >& x0,
  std::vector<Contact> & contact)
{
  TrigMesh * surf = trigm[m0];
  ElementMesh * em = em_[m0];
  ElementMesh * em1 = em_[m1];
  double dx = em->X[1][2] - em->X[0][2];
  float eps = 0;
  Eigen::Vector3d trig[3], trig0[3];
  Eigen::Vector3d trigNormal;

  for (int l = 0; l < 3; l++){
    int vi = surf->vidx[surf->t[t0][l]];
    trig[l] = em->x[vi];
    trig0[l] = x0[m0][vi];
  }
  trigNormal = (trig[1] - trig[0]).cross(trig[2] - trig[0]);
  trigNormal.normalize();
  for (int l = 0; l < 3; l++){
    trig0[l] -= trigNormal * eps;
    trig[l] += trigNormal * eps;
  }
  Eigen::Vector3d x_1 = em1->x[v1];
  Eigen::Vector3d x_0 = x0[m1][v1];

  InterResult result = PtTrigIntersect(x_0, x_1, trig0, trig);
  if (result.intersect){
    Contact c;
    c.m[0] = m1;
    c.m[1] = m0;
    c.v1 = v1;
    c.normal = result.n;
    c.x = (1 - result.t) * x_0 + result.t * x_1;

    //quad vertex indices in the finite element mesh.
    for (int l = 0; l < 4; l++){
      int vidx = em->e[surf->ei[t0]]->at(HexFaces[surf->fi[t0]][l]);
      c.v2[l] = vidx;
    }

    //corner of quad.
    Eigen::Vector3d x0 = em->x[c.v2[0]];
    Eigen::Vector3d x1 = em->x[c.v2[1]];
    Eigen::Vector3d x2 = em->x[c.v2[2]];

    Eigen::Vector3d a1 = (x1 - x0);
    double l1 = a1.norm();
    a1 /= l1;
    Eigen::Vector3d a2 = (x2 - x0);
    double l2 = a2.norm();
    Eigen::Vector3d a3 = a1.cross(a2).normalized();
    a2 = a3.cross(a1);
    Eigen::Vector3d disp = c.x - x0;
    c.alpha[0] = (float)(a1.dot(disp) / l1);
    c.alpha[1] = (float)(a2.dot(disp) / l2);
    contact.push_back(c);
    return true;
  }
  return false;
}

bool addContact(ElementMesh * em,
  TrigMesh * trigm,
  int m0, int t0, int v1,
  const std::vector<Eigen::Vector3d> & x0,
  std::vector<Contact> & contact)
{
  TrigMesh * surf = trigm;
  double dx = em->X[1][2] - em->X[0][2];
  float eps = 0;
  Eigen::Vector3d trig[3], trig0[3];
  Eigen::Vector3d trigNormal;

  for (int l = 0; l < 3; l++){
    int vi = surf->vidx[surf->t[t0][l]];
    trig[l] = em->x[vi];
    trig0[l] = x0[vi];
  }
  trigNormal = (trig[1] - trig[0]).cross(trig[2] - trig[0]);
  trigNormal.normalize();
  for (int l = 0; l < 3; l++){
    trig0[l] -= trigNormal * eps;
    trig[l] += trigNormal * eps;
  }
  Eigen::Vector3d x_1 = em->x[v1];
  Eigen::Vector3d x_0 = x0[v1];

  InterResult result = PtTrigIntersect(x_0, x_1, trig0, trig);
  if (result.intersect){
    Contact c;
    c.m[0] = m0;
    c.m[1] = m0;
    c.v1 = v1;
    c.normal = result.n;
    c.x = (1 - result.t) * x_0 + result.t * x_1;

    //quad vertex indices in the finite element mesh.
    for (int l = 0; l < 4; l++){
      int vidx = em->e[surf->ei[t0]]->at(HexFaces[surf->fi[t0]][l]);
      c.v2[l] = vidx;
    }

    //corner of quad.
    Eigen::Vector3d x0 = em->x[c.v2[0]];
    Eigen::Vector3d x1 = em->x[c.v2[1]];
    Eigen::Vector3d x2 = em->x[c.v2[2]];

    Eigen::Vector3d a1 = (x1 - x0);
    double l1 = a1.norm();
    a1 /= l1;
    Eigen::Vector3d a2 = (x2 - x0);
    double l2 = a2.norm();
    Eigen::Vector3d a3 = a1.cross(a2).normalized();
    a2 = a3.cross(a1);
    Eigen::Vector3d disp = c.x - x0;
    c.alpha[0] = (float)(a1.dot(disp) / l1);
    c.alpha[1] = (float)(a2.dot(disp) / l2);
    contact.push_back(c);
    return true;
  }
  return false;
}

bool detectCollision(std::vector<ElementMesh *> & em_,
  std::vector<TrigMesh *> trigm,
  const std::vector<std::vector<Eigen::Vector3d> >& x0,
  std::vector<Contact> & contact)
{
  contact.clear();
  for (size_t mi = 0; mi < em_.size(); mi++){
    updateSurfVert(em_[mi]->x, trigm[mi], trigm[mi]->vidx);
  }

  //detect bound collision
  for (size_t mi = 0; mi < em_.size(); mi++){
    for (int j = 0; j < em_[mi]->x.size(); j++){
      if (em_[mi]->fixedDof[3 * j]){
        continue;
      }
      if (em_[mi]->x[j][1] < em_[mi]->lb[j][1]){
        Contact c;
        c.m[0] = (int)mi;
        c.m[1] = -1;
        c.v1 = j;
        c.normal = Eigen::Vector3d(0, 1, 0);
        c.x = em_[mi]->x[j];
        contact.push_back(c);
      }
    }
  }

  //detect vertex-quad collision
  if (em_.size() >= 2){
    int midx0 = 0;
    int midx1 = 1;
    TrigMesh * surf = trigm[midx0];
    TrigMesh * surf1 = trigm[midx1];
    ElementMesh * m = em_[midx0];
    ElementMesh * m1 = em_[midx1];
    for (int j = 0; j < surf1->v.size(); j++){
      int vidx = surf1->vidx[j];
      if (m1->fixedDof[3 * vidx]){
        continue;
      }
      for (int k = 0; k < surf->t.size(); k++){
        bool hasContact = addContact(em_, trigm, midx0, midx1, k, vidx, x0, contact);
        if (hasContact){
          break;
        }
      }
    }

    midx0 = 1;
    midx1 = 0;
    surf = trigm[midx0];
    surf1 = trigm[midx1];
    m = em_[midx0];
    m1 = em_[midx1];
    for (int j = 0; j < surf1->v.size(); j++){
      int vidx = surf1->vidx[j];
      if (vidx < 0 || m1->fixedDof[3 * j]){
        continue;
      }
      for (int k = 0; k < surf->t.size(); k++){
        bool hasContact = addContact(em_, trigm, midx0, midx1, k, j, x0, contact);
        if (hasContact){
          break;
        }
      }
    }

  }
  return contact.size()>0;
}
